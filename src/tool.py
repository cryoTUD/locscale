"""
Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the
program 'LocScale2' written by the Author(s).

Copyright (C) 2026 Alok Bharadwaj and Arjen J. Jakobi

GUI for the LocScale2 feature-enhance pipeline. Structure follows LocScale-SURFER: a
ToolInstance holding a MainToolWindow, volumes chosen with ModelMenuButton, and the pipeline
on a QThread so ChimeraX stays responsive.
"""
import json
import os
import re

import numpy as np

from chimerax.core.models import Model, Surface
from chimerax.core.tools import ToolInstance
from chimerax.map import Volume
from chimerax.ui import MainToolWindow
from chimerax.ui.widgets import CollapsiblePanel, ModelMenuButton, vertical_layout
from Qt.QtCore import QThread, Signal
from Qt.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QFrame, QGroupBox,
                          QHBoxLayout, QHeaderView, QLabel, QLineEdit, QProgressBar,
                          QPushButton, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout)


# Give every boxed section (QGroupBox) and collapsible section a light rounded border with a
# bold title, so the panels read as distinct parts. palette(mid) tracks the light/dark theme.
_SECTION_STYLE = """
QGroupBox, CollapsiblePanel {
    border: 1px solid palette(mid);
    border-radius: 6px;
    margin-top: 8px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 8px;
    padding: 0 4px;
    font-weight: bold;
}
"""


from .utils import (
    round_up_proper,
    round_up_to_even,
)

with open(os.path.join(os.path.dirname(__file__), "data", "help_info.json")) as _f:
    help_info = json.load(_f)


def _is_hex_color(text):
    """True for a '#rrggbb' string."""
    return bool(re.fullmatch(r"#[0-9a-fA-F]{6}", (text or "").strip()))


def _reverse_palette_spec(spec):
    """Reverse a 'value,color:value,color:...' palette's colour direction.

    Keeps the (ascending) pVDDT values in place and reverses only the colour order, so e.g.
    blue-green-red becomes red-green-blue.
    """
    pairs = [p.split(",", 1) for p in spec.split(":")]
    values = [v for v, _ in pairs]
    colors = [c for _, c in pairs]
    return ":".join("{},{}".format(v, c) for v, c in zip(values, reversed(colors)))


class PipelineWorker(QThread):
    """Runs the pipeline off the UI thread.

    Nothing here may touch ChimeraX models: results come back through `completed` and are
    turned into volumes on the UI thread. The signal is deliberately not called `finished`,
    which is QThread's own.
    """

    status = Signal(str)
    progress = Signal(str, int, int)
    completed = Signal(object)
    cancelled = Signal()
    failed = Signal(str)

    def __init__(self, kwargs, run_kind="feature_enhance"):
        super().__init__()
        self._kwargs = kwargs
        self._run_kind = run_kind
        self._cancel_requested = False

    def cancel(self):
        """Ask the run to stop. Safe to call from the UI thread: it only sets a flag."""
        self._cancel_requested = True

    def run(self):
        from .pipeline import Cancelled, run_amplitude_scaling, run_feature_enhance

        run_fn = (run_amplitude_scaling if self._run_kind == "amplitude_scaling"
                  else run_feature_enhance)

        # The pipeline's callbacks are the only points at which it yields to us, so they are
        # also where cancellation is honoured.
        def status(message):
            if self._cancel_requested:
                raise Cancelled()
            self.status.emit(message)

        def progress(stage, done, total):
            if self._cancel_requested:
                raise Cancelled()
            self.progress.emit(stage, done, total)

        try:
            results = run_fn(
                status_callback=status, progress_callback=progress, **self._kwargs)
        except Cancelled:
            self.cancelled.emit()
            return
        except Exception as exc:
            import traceback
            self.failed.emit("{}: {}\n\n{}".format(type(exc).__name__, exc, traceback.format_exc()))
            return
        self.completed.emit(results)


class LocScale2Tool(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "LocScale-FEM"
        self._worker = None
        self._template = None
        self._result_kind = "feature_enhance"   # which display the running worker feeds
        self._verify_reference = None     # set by 'locscale2 verify' for the CC report
        self._noise_box_model = None      # parent Model holding the box surfaces
        self._noise_edited = False        # user has hand-edited the boxes/window
        self._populating = False          # guard while filling the table programmatically

        self.tool_window = MainToolWindow(self)
        parent = self.tool_window.ui_area
        parent.setStyleSheet(_SECTION_STYLE)
        layout = vertical_layout(parent, margins=(5, 5, 5, 5), spacing=6)
        layout.addWidget(self._inputs_panel(parent))
        layout.addWidget(self._symmetry_panel(parent))
        layout.addWidget(self._options_panel(parent))
        layout.addWidget(self._run_panel(parent))
        layout.addWidget(self._pvddt_color_panel(parent))
        layout.addWidget(self._citation_label(parent))
        layout.addStretch(1)
        self.tool_window.manage("side")
        self._populate_noise_defaults()

        # # Show the options expanded. Done after manage() so the content area's sizeHint is
        # # settled; setChecked keeps the disclosure button in step, otherwise the first
        # # click would try to expand an already-expanded panel.
        # self._options.toggle_button.setChecked(True)
        # self._options.toggle_panel_display(True)

    # ------------------------------------------------------------ panels

    def _row(self, parent, label_text, tooltip, widget):
        row = QFrame(parent)
        hl = QHBoxLayout(row)
        hl.setContentsMargins(0, 0, 0, 0)
        label = QLabel(label_text, row)
        label.setToolTip(tooltip)
        hl.addWidget(label)
        hl.addWidget(widget)
        hl.addStretch(1)
        return row

    def _group(self, parent, title):
        """A titled, bordered box with a vertical layout; returns (groupbox, layout)."""
        gb = QGroupBox(title, parent)
        lay = QVBoxLayout(gb)
        lay.setContentsMargins(8, 4, 8, 8)
        return gb, lay

    def _inputs_panel(self, parent):
        frame, panel = self._group(parent, "Inputs")

        # Input unsharpened map
        self._map_menu = ModelMenuButton(self.session, class_filter=Volume)
        self._map_menu.value_changed.connect(self._on_input_map_changed)
        panel.addWidget(self._row(frame, "Input map:", help_info["input_map_help"], self._map_menu))

        # Alternatively, two half maps that are averaged into the working map.
        self._half1_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        self._half1_menu.value_changed.connect(self._on_input_map_changed)
        panel.addWidget(self._row(frame, "Half map 1 (optional):",
                                  "First half map; averaged with half map 2 to form the "
                                  "working map. Overrides 'Input map' when both are set.",
                                  self._half1_menu))
        self._half2_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        self._half2_menu.value_changed.connect(self._on_input_map_changed)
        panel.addWidget(self._row(frame, "Half map 2 (optional):",
                                  "Second half map; averaged with half map 1 to form the "
                                  "working map.", self._half2_menu))
        half_note = QLabel("<i>Give an input map, or two half maps to be averaged.</i>", frame)
        half_note.setWordWrap(True)
        panel.addWidget(half_note)

        # Optional mask
        self._mask_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        panel.addWidget(self._row(frame, "Mask (optional):", help_info["mask_help"], self._mask_menu))
        note = QLabel("<i>If no mask is given, an FDR mask is computed and returned.</i>", frame)
        note.setWordWrap(True)
        panel.addWidget(note)

        # FDR noise-box controls (collapsible) live inside Inputs.
        panel.addWidget(self._noise_panel(frame))
        return frame

    def _active_input_volume(self):
        """The Volume whose grid drives noise-box defaults/surfaces: the input map, else the
        first half map when a half-map pair is chosen."""
        if self._map_menu.value is not None:
            return self._map_menu.value
        if self._half1_menu.value is not None and self._half2_menu.value is not None:
            return self._half1_menu.value
        return None

    def _resolve_working_map(self):
        """Return (emmap, apix, template_volume), or None after setting an error message.

        Two half maps (both set) are averaged into the working map; otherwise the single
        input map is used.
        """
        input_map = self._map_menu.value
        half1, half2 = self._half1_menu.value, self._half2_menu.value
        if half1 is not None or half2 is not None:
            if half1 is None or half2 is None:
                self._status_label.setText(
                    "<font color='red'>Select both half maps, or neither.</font>")
                return None
            if half1 is half2:
                self._status_label.setText(
                    "<font color='red'>The two half maps are the same volume.</font>")
                return None
            h1 = np.asarray(half1.data.full_matrix(), dtype=np.float32)
            h2 = np.asarray(half2.data.full_matrix(), dtype=np.float32)
            if h1.shape != h2.shape:
                self._status_label.setText(
                    "<font color='red'>Half map shapes {} and {} differ.</font>".format(
                        h1.shape, h2.shape))
                return None
            return 0.5 * (h1 + h2), float(half1.data.step[0]), half1
        if input_map is None:
            self._status_label.setText(
                "<font color='red'>Select an input map, or two half maps.</font>")
            return None
        return input_map.data.full_matrix(), float(input_map.data.step[0]), input_map

    def _symmetry_panel(self, parent):
        """Collapsible 'Symmetrise output?' section: point-group and helical symmetry."""
        self._symmetry_section = panel = CollapsiblePanel(parent, title="Symmetrise output?")
        frame = panel.content_area
        box = frame.layout()

        self._point_group_symmetry_menu = QLineEdit(frame)
        box.addWidget(self._row(frame, "Point group symmetry:",
                                help_info["point_group_symmetry_help"],
                                self._point_group_symmetry_menu))

        self._helical_symmetry_checkbox = QCheckBox("Helical symmetry", frame)
        self._helical_symmetry_checkbox.setToolTip(help_info["helical_symmetry_help"])
        self._helical_symmetry_checkbox.toggled.connect(self._helical_symmetry_toggled)
        self._helical_symmetry_checkbox.setChecked(False)
        box.addWidget(self._row(frame, "Helical symmetry:", help_info["helical_symmetry_help"],
                                self._helical_symmetry_checkbox))

        # Helical rise/twist inputs, shown only while the helical checkbox is ticked.
        self._helical_symmetry_panel = QFrame(frame)
        helical_layout = QVBoxLayout(self._helical_symmetry_panel)
        helical_layout.setContentsMargins(0, 0, 0, 0)
        self._helical_rise_menu = QLineEdit(self._helical_symmetry_panel, placeholderText="0.0")
        helical_layout.addWidget(self._row(self._helical_symmetry_panel, "Rise (A):",
                                           help_info["helical_symmetry_help"], self._helical_rise_menu))
        self._helical_twist_menu = QLineEdit(self._helical_symmetry_panel, placeholderText="0.0")
        helical_layout.addWidget(self._row(self._helical_symmetry_panel, "Twist (deg):",
                                           help_info["helical_symmetry_help"], self._helical_twist_menu))
        self._helical_symmetry_panel.setVisible(False)
        box.addWidget(self._helical_symmetry_panel)
        return panel

    def _helical_symmetry_toggled(self, checked):
        """Show the rise/twist inputs only while helical symmetry is requested."""
        self._helical_symmetry_panel.setVisible(checked)
        # Re-fit the collapsible so the newly shown rows are not clipped.
        if self._symmetry_section.shown:
            self._symmetry_section.resize_panel(True)

    def _noise_panel(self, parent):
        panel = CollapsiblePanel(parent, title="FDR mask - noise boxes")
        frame = panel.content_area
        box = frame.layout()

        note = QLabel("<i>Used only when no mask is supplied. Boxes are (x, y, z) voxel "
                      "centres; their pooled voxels give the noise mean/variance.</i>", frame)
        note.setWordWrap(True)
        box.addWidget(note)

        self._noise_table = QTableWidget(0, 3, frame)
        self._noise_table.setHorizontalHeaderLabels(["x", "y", "z"])
        self._noise_table.verticalHeader().setVisible(False)
        self._noise_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._noise_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._noise_table.setMaximumHeight(160)
        self._noise_table.itemChanged.connect(self._on_noise_edited)
        box.addWidget(self._noise_table)

        buttons = QFrame(frame)
        bl = QHBoxLayout(buttons); bl.setContentsMargins(0, 0, 0, 0)
        add_btn = QPushButton("Add", buttons); add_btn.clicked.connect(self._add_noise_row)
        rm_btn = QPushButton("Remove", buttons); rm_btn.clicked.connect(self._remove_noise_row)
        reset_btn = QPushButton("Reset", buttons)
        reset_btn.clicked.connect(lambda: self._populate_noise_defaults(force=True))
        bl.addWidget(add_btn); bl.addWidget(rm_btn); bl.addWidget(reset_btn); bl.addStretch(1)
        box.addWidget(buttons)

        self._noise_window_spin = QSpinBox(frame)
        self._noise_window_spin.setRange(4, 512)
        self._noise_window_spin.setValue(20)
        self._noise_window_spin.valueChanged.connect(self._on_noise_window_changed)
        box.addWidget(self._row(frame, "Noise box window (px):",
                                "Cube edge, in pixels, sampled around each box centre. "
                                "Default: 10% of the map edge or 20 px, whichever is larger.",
                                self._noise_window_spin))

        self._noise_show_check = QCheckBox("Visualize noise boxes", frame)
        self._noise_show_check.toggled.connect(self._toggle_noise_surfaces)
        box.addWidget(self._noise_show_check)
        return panel

    # ------------------------------------------------------------ noise-box logic

    def _add_noise_row(self):
        self._noise_edited = True
        self._append_noise_row(0.0, 0.0, 0.0)
        self._refresh_noise_surfaces()

    def _remove_noise_row(self):
        rows = sorted({i.row() for i in self._noise_table.selectedIndexes()}, reverse=True)
        if not rows:
            rows = [self._noise_table.rowCount() - 1]
        self._noise_edited = True
        for r in rows:
            if r >= 0:
                self._noise_table.removeRow(r)
        self._refresh_noise_surfaces()

    def _append_noise_row(self, x, y, z):
        r = self._noise_table.rowCount()
        self._noise_table.insertRow(r)
        for c, val in enumerate((x, y, z)):
            self._noise_table.setItem(r, c, QTableWidgetItem("{:.1f}".format(val)))

    def _read_noise_boxes(self):
        """Valid (x, y, z) rows from the table; incomplete/invalid rows are skipped."""
        boxes = []
        for r in range(self._noise_table.rowCount()):
            try:
                vals = [float(self._noise_table.item(r, c).text()) for c in range(3)]
            except (AttributeError, ValueError):
                continue
            boxes.append(tuple(vals))
        return boxes

    def _populate_noise_defaults(self, force=False):
        """Fill the window size and the four edge-patch boxes for the current input map."""
        if self._noise_edited and not force:
            return
        volume = self._active_input_volume()
        if volume is None:
            return
        from .pipeline import default_noise_window_size
        from .include.mapops import default_noise_box_centers
        xs, ys, zs = volume.data.size          # grid size (x, y, z); no data load
        shape = (zs, ys, xs)                    # match full_matrix() array order
        window = default_noise_window_size(shape)
        self._populating = True
        try:
            self._noise_window_spin.setValue(window)
            self._noise_table.setRowCount(0)
            for x, y, z in default_noise_box_centers(shape, window):
                self._append_noise_row(x, y, z)
        finally:
            self._populating = False
        self._noise_edited = False
        self._refresh_noise_surfaces()

    def _on_input_map_changed(self):
        if getattr(self, "_noise_table", None) is None:
            return                            # signal fired before the panel was built
        self._populate_noise_defaults()
        if self._noise_edited:
            self._refresh_noise_surfaces()   # re-place existing boxes on the new map

    def _on_noise_edited(self, *args):
        if self._populating:
            return
        self._noise_edited = True
        self._refresh_noise_surfaces()

    def _on_noise_window_changed(self, *args):
        if not self._populating:
            self._noise_edited = True
        self._refresh_noise_surfaces()

    # ------------------------------------------------------------ noise-box surfaces

    def _toggle_noise_surfaces(self, checked):
        self._refresh_noise_surfaces()

    def _clear_noise_surfaces(self):
        if self._noise_box_model is not None and not self._noise_box_model.deleted:
            self.session.models.close([self._noise_box_model])
        self._noise_box_model = None

    def _refresh_noise_surfaces(self):
        self._clear_noise_surfaces()
        if not self._noise_show_check.isChecked():
            return
        volume = self._active_input_volume()
        boxes = self._read_noise_boxes()
        if volume is None or not boxes:
            return
        window = self._noise_window_spin.value()
        parent = Model("Noise boxes", self.session)
        for i, center in enumerate(boxes):
            parent.add([self._build_box_surface(center, window, volume, i)])
        volume.add([parent])                 # nest under the map so it inherits its placement
        self._noise_box_model = parent

    def _build_box_surface(self, center, window, volume, index):
        h = 0.5 * window
        cx, cy, cz = center
        signs = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                 (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
        ijk = [(cx + sx * h, cy + sy * h, cz + sz * h) for sx, sy, sz in signs]
        vertices = np.array([volume.data.ijk_to_xyz(p) for p in ijk], dtype=np.float32)
        triangles = np.array([
            (0, 1, 2), (0, 2, 3),            # -z
            (4, 6, 5), (4, 7, 6),            # +z
            (0, 4, 5), (0, 5, 1),            # -y
            (1, 5, 6), (1, 6, 2),            # +x
            (2, 6, 7), (2, 7, 3),            # +y
            (3, 7, 4), (3, 4, 0),            # -x
        ], dtype=np.int32)
        from chimerax.surface import calculate_vertex_normals
        normals = calculate_vertex_normals(vertices, triangles)
        surface = Surface("box {}".format(index + 1), self.session)
        surface.set_geometry(vertices, normals, triangles)
        surface.color = np.array([255, 210, 40, 120], dtype=np.uint8)   # translucent amber
        return surface

    # ------------------------------------------------------------ pVDDT colouring

    def _pvddt_color_panel(self, parent):
        from .pipeline import PVDDT_PALETTES, PVDDT_CUSTOM_STOPS
        self._pvddt_panel = panel = CollapsiblePanel(parent, title="Colour by pVDDT")
        frame = panel.content_area
        box = frame.layout()

        note = QLabel("<i>Colours a displayed map's surface by a pVDDT map -- no feature "
                      "enhancement is run. The map to colour must be shown as a surface.</i>",
                      frame)
        note.setWordWrap(True)
        box.addWidget(note)

        self._pvddt_map_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        box.addWidget(self._row(frame, "Map to colour:", "The map whose surface is coloured.",
                                self._pvddt_map_menu))

        self._pvddt_value_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        box.addWidget(self._row(frame, "pVDDT map:", "The pVDDT map providing the colours.",
                                self._pvddt_value_menu))

        self._pvddt_scheme_combo = QComboBox(frame)
        self._pvddt_scheme_combo.addItems(list(PVDDT_PALETTES.keys()))
        self._pvddt_scheme_combo.addItem("Custom...")
        self._pvddt_scheme_combo.currentIndexChanged.connect(self._on_pvddt_scheme_changed)
        box.addWidget(self._row(frame, "Colour scheme:",
                                "Diverging palette over pVDDT -100..+100.",
                                self._pvddt_scheme_combo))

        self._pvddt_reverse_check = QCheckBox("Reverse colour direction", frame)
        self._pvddt_reverse_check.setToolTip(
            "Flip the palette, e.g. blue-green-red becomes red-green-blue.")
        box.addWidget(self._pvddt_reverse_check)

        # Custom-palette editor, revealed only when "Custom..." is chosen. Kept as a plain
        # frame (not a nested CollapsiblePanel): a collapsible inside a collapsible freezes
        # the outer panel's height and clips these rows.
        self._pvddt_custom_frame = QFrame(frame)
        cbox = QVBoxLayout(self._pvddt_custom_frame)
        cbox.setContentsMargins(0, 0, 0, 0)
        cbox.addWidget(QLabel("<b>Custom palette</b>", self._pvddt_custom_frame))
        self._pvddt_custom_edits = []
        for value, hexcolor in PVDDT_CUSTOM_STOPS:
            edit = QLineEdit(hexcolor, self._pvddt_custom_frame)
            self._pvddt_custom_edits.append((value, edit))
            cbox.addWidget(self._row(self._pvddt_custom_frame, "pVDDT {:+d}:".format(value),
                                     "Hex colour (#rrggbb) for pVDDT {:+d}.".format(value),
                                     edit))
        self._pvddt_custom_frame.setVisible(False)
        box.addWidget(self._pvddt_custom_frame)

        self._pvddt_apply_button = QPushButton("Apply colour", frame)
        self._pvddt_apply_button.clicked.connect(self._apply_pvddt_color)
        box.addWidget(self._pvddt_apply_button)
        return panel

    def _on_pvddt_scheme_changed(self, index):
        is_custom = self._pvddt_scheme_combo.currentText() == "Custom..."
        self._pvddt_custom_frame.setVisible(is_custom)
        # The outer panel froze its height when it was expanded; re-fit so the revealed rows
        # are not clipped (same fix as _gpu_toggled).
        if self._pvddt_panel.shown:
            self._pvddt_panel.resize_panel(True)

    def _pvddt_palette_spec(self):
        """Palette spec for the current scheme, or None if a custom hex is invalid."""
        from .pipeline import PVDDT_PALETTES
        name = self._pvddt_scheme_combo.currentText()
        if name != "Custom...":
            spec = PVDDT_PALETTES[name]
        else:
            stops = []
            for value, edit in self._pvddt_custom_edits:
                hexcolor = edit.text().strip()
                if not _is_hex_color(hexcolor):
                    return None
                stops.append("{},{}".format(value, hexcolor))
            spec = ":".join(stops)
        if self._pvddt_reverse_check.isChecked():
            spec = _reverse_palette_spec(spec)
        return spec

    def _apply_pvddt_color(self):
        map_vol = self._pvddt_map_menu.value
        pvddt_vol = self._pvddt_value_menu.value
        if map_vol is None or pvddt_vol is None:
            self._status_label.setText(
                "<font color='red'>Select a map to colour and a pVDDT map.</font>")
            return
        spec = self._pvddt_palette_spec()
        if spec is None:
            self._status_label.setText(
                "<font color='red'>Custom palette has an invalid hex colour (use #rrggbb).</font>")
            return
        from chimerax.core.commands import run
        try:
            run(self.session, "color sample #{} map #{} palette {}".format(
                map_vol.id_string, pvddt_vol.id_string, spec))
        except Exception as exc:
            self._status_label.setText("<font color='red'>Colouring failed; see the log.</font>")
            self.session.logger.error("LocScale2 colour-by-pVDDT failed: {}".format(exc))
            return
        self._status_label.setText("Coloured {} by pVDDT.".format(map_vol.name))

    def _options_panel(self, parent):
        self._options = panel = CollapsiblePanel(parent, title="Advanced Options")
        frame = panel.content_area
        # CollapsiblePanel already installs a QVBoxLayout on content_area; reuse it.
        options = frame.layout()

        # --- EMmerNet (prediction) options -------------------------------------
        emm, emm_l = self._group(frame, "EMmerNet options")
        self._model_combo = QComboBox(emm)
        from .emmernet import available_models
        self._model_combo.addItems(available_models())
        emm_l.addWidget(self._row(emm, "Model:", help_info["model_help"], self._model_combo))

        self._mc_spin = QSpinBox(emm); self._mc_spin.setRange(1, 100); self._mc_spin.setValue(15)
        emm_l.addWidget(self._row(emm, "Monte-Carlo iterations:",
                                  help_info["monte_carlo_help"], self._mc_spin))

        self._batch_spin = QSpinBox(emm); self._batch_spin.setRange(1, 64); self._batch_spin.setValue(8)
        emm_l.addWidget(self._row(emm, "Batch size:", help_info["batch_size_help"], self._batch_spin))

        # Cube size is deliberately not exposed: EMmerNet's padding and crops assume 32.
        self._stride_spin = QSpinBox(emm); self._stride_spin.setRange(1, 32); self._stride_spin.setValue(16)
        emm_l.addWidget(self._row(emm, "Cube stride:", help_info["stride_help"], self._stride_spin))

        self._gpu_check = QCheckBox("Use GPU when available", emm)
        self._gpu_check.setChecked(True)
        self._gpu_check.setToolTip(help_info["gpu_help"])
        self._gpu_check.toggled.connect(self._gpu_toggled)
        emm_l.addWidget(self._gpu_check)

        from .emmernet import available_gpus
        self._gpus = available_gpus()
        self._gpu_combo = QComboBox(emm)
        for index, name in self._gpus:
            self._gpu_combo.addItem("{}: {}".format(index, name), index)
        if not self._gpus:
            # MPS and CPU have no device index to pick, so leave the control in place but
            # inert rather than implying a choice that does not exist.
            self._gpu_combo.addItem("Using CPU", None)
            self._gpu_combo.setEnabled(False)
        self._gpu_row = self._row(emm, "GPU:", help_info["gpu_id_help"], self._gpu_combo)
        emm_l.addWidget(self._gpu_row)

        # --- Amplitude-scaling options -----------------------------------------
        scl, scl_l = self._group(frame, "Scaling options")
        self._window_spin = QSpinBox(scl); self._window_spin.setRange(11, 45); self._window_spin.setValue(25)
        scl_l.addWidget(self._row(scl, "Scaling window:", help_info["window_help"], self._window_spin))

        self._chunk_spin = QSpinBox(scl)
        self._chunk_spin.setRange(128, 65536); self._chunk_spin.setSingleStep(512)
        self._chunk_spin.setValue(4096)
        scl_l.addWidget(self._row(scl, "Scaling chunk:", help_info["chunk_help"], self._chunk_spin))
        scl_l.addStretch(1)      # keep the shorter box top-aligned beside EMmerNet options

        # EMmerNet and Scaling options sit side by side (two-column grid).
        top_row = QHBoxLayout()
        top_row.setContentsMargins(0, 0, 0, 0)
        top_row.addWidget(emm)
        top_row.addWidget(scl)
        options.addLayout(top_row)

        # --- Reference-based (model-based) LocScale ----------------------------
        ref, ref_l = self._group(frame, "Run model-based LocScale")
        ref_note = QLabel("<i>Scale the input map to a reference map only -- skips feature "
                          "enhancement and pVDDT.</i>", ref)
        ref_note.setWordWrap(True)
        ref_l.addWidget(ref_note)

        ref_row = QFrame(ref)
        ref_layout = QHBoxLayout(ref_row)
        ref_layout.setContentsMargins(0, 0, 0, 0)
        self._reference_menu = ModelMenuButton(
            self.session, class_filter=Volume,
            no_value_button_text="No model chosen", no_value_menu_text="None",
            autoselect="none")
        self._run_locscale_button = QPushButton("Run LocScale", ref_row)
        self._run_locscale_button.setToolTip(
            "Amplitude-scale the input map to the reference map only -- no feature "
            "enhancement, no pVDDT.")
        self._run_locscale_button.clicked.connect(self._run_locscale_clicked)
        ref_layout.addWidget(QLabel("Reference map:", ref_row))
        ref_layout.addWidget(self._reference_menu)
        ref_layout.addWidget(self._run_locscale_button)
        ref_layout.addStretch(1)
        ref_l.addWidget(ref_row)
        options.addWidget(ref)
        return panel

    _CITATION_DOI = "10.1038/s41467-026-75327-8"

    def _citation_label(self, parent):
        label = QLabel(
            "<i>If LocScale-FEM is useful in your work, please cite: "
            "<a href='https://doi.org/{doi}'>doi:{doi}</a></i>".format(doi=self._CITATION_DOI),
            parent)
        label.setWordWrap(True)
        label.setOpenExternalLinks(True)
        return label

    def _run_panel(self, parent):
        frame = QFrame(parent)
        panel = vertical_layout(frame, margins=(0, 8, 0, 0))

        buttons = QFrame(frame)
        button_row = QHBoxLayout(buttons)
        button_row.setContentsMargins(0, 0, 0, 0)

        self._run_button = QPushButton("Run feature enhancement", buttons)
        self._run_button.clicked.connect(self._run_clicked)
        button_row.addWidget(self._run_button)

        self._cancel_button = QPushButton("Cancel", buttons)
        self._cancel_button.setEnabled(False)
        self._cancel_button.setToolTip(help_info["cancel_help"])
        self._cancel_button.clicked.connect(self._cancel_clicked)
        button_row.addWidget(self._cancel_button)

        panel.addWidget(buttons)

        self._progress_bar = QProgressBar(frame)
        self._progress_bar.setVisible(False)
        panel.addWidget(self._progress_bar)

        self._status_label = QLabel("", frame)
        self._status_label.setWordWrap(True)
        panel.addWidget(self._status_label)
        return frame

    # ------------------------------------------------------------ actions

    def _run_clicked(self):
        resolved = self._resolve_working_map()
        if resolved is None:
            return
        emmap, apix, template = resolved
        mask_volume = self._mask_menu.value
        if mask_volume is template:
            self._status_label.setText("<font color='red'>Map and mask are the same volume.</font>")
            return

        window_size_angstroms = self._window_spin.value()
        window_size_pix = round_up_to_even(window_size_angstroms / apix)

        mask = mask_volume.data.full_matrix() if mask_volume is not None else None
        if mask is not None and mask.shape != emmap.shape:
            self._status_label.setText(
                "<font color='red'>Mask shape {} does not match map {}.</font>".format(
                    mask.shape, emmap.shape))
            return

        self._template = template
        self._result_kind = "feature_enhance"
        self._set_running(True)

        # Get the point group symmetry and helical symmetry parameters.
        # An empty field means "no symmetry": placeholderText is only a hint, not a value,
        # so .text() is "" when the user leaves it blank -- default that to C1.
        pg = self._point_group_symmetry_menu.text().strip() or "C1"
        helical_symmetry = self._helical_symmetry_checkbox.isChecked()
        rise = float(self._helical_rise_menu.text()) if helical_symmetry else None
        twist = float(self._helical_twist_menu.text()) if helical_symmetry else None

        noise_boxes = self._read_noise_boxes() or None
        self._worker = PipelineWorker(dict(
            emmap=emmap,
            apix=apix,
            mask=mask,
            noise_boxes=noise_boxes,                       # used only when mask is None
            noise_window_size=self._noise_window_spin.value(),
            model_type=self._model_combo.currentText(),
            monte_carlo_iterations=self._mc_spin.value(),
            batch_size=self._batch_spin.value(),
            stride=self._stride_spin.value(),
            window_size=window_size_pix,
            scaling_chunk=self._chunk_spin.value(),
            use_gpu=self._gpu_check.isChecked(),
            gpu_id=self._selected_gpu_id(),
            pg=pg,
            rise=rise,
            twist=twist
        ))

        self._worker.status.connect(self._on_status)
        self._worker.progress.connect(self._on_progress)
        self._worker.completed.connect(self._on_completed)
        self._worker.cancelled.connect(self._on_cancelled)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _run_locscale_clicked(self):
        """Amplitude scaling only, against the chosen reference map (skips EMmerNet)."""
        reference_map = self._reference_menu.value
        if reference_map is None:
            self._status_label.setText(
                "<font color='red'>Select a reference map for LocScale.</font>")
            return

        resolved = self._resolve_working_map()
        if resolved is None:
            return
        emmap, apix, template = resolved
        mask_volume = self._mask_menu.value

        if reference_map is template:
            self._status_label.setText(
                "<font color='red'>Reference and input map are the same volume.</font>")
            return
        if mask_volume is template:
            self._status_label.setText("<font color='red'>Map and mask are the same volume.</font>")
            return

        reference = reference_map.data.full_matrix()
        if reference.shape != emmap.shape:
            self._status_label.setText(
                "<font color='red'>Reference shape {} does not match map {}.</font>".format(
                    reference.shape, emmap.shape))
            return

        window_size_angstroms = self._window_spin.value()
        window_size_pix = round_up_to_even(window_size_angstroms / apix)

        mask = mask_volume.data.full_matrix() if mask_volume is not None else None
        if mask is not None and mask.shape != emmap.shape:
            self._status_label.setText(
                "<font color='red'>Mask shape {} does not match map {}.</font>".format(
                    mask.shape, emmap.shape))
            return

        self._template = template
        self._result_kind = "amplitude_scaling"
        self._set_running(True)

        pg = self._point_group_symmetry_menu.text().strip() or "C1"
        helical_symmetry = self._helical_symmetry_checkbox.isChecked()
        rise = float(self._helical_rise_menu.text()) if helical_symmetry else None
        twist = float(self._helical_twist_menu.text()) if helical_symmetry else None

        noise_boxes = self._read_noise_boxes() or None
        self._worker = PipelineWorker(dict(
            emmap=emmap,
            apix=apix,
            reference=reference,
            mask=mask,
            noise_boxes=noise_boxes,                       # used only when mask is None
            noise_window_size=self._noise_window_spin.value(),
            window_size=window_size_pix,
            scaling_chunk=self._chunk_spin.value(),
            use_gpu=self._gpu_check.isChecked(),
            gpu_id=self._selected_gpu_id(),
            pg=pg,
            rise=rise,
            twist=twist,
        ), run_kind="amplitude_scaling")

        self._worker.status.connect(self._on_status)
        self._worker.progress.connect(self._on_progress)
        self._worker.completed.connect(self._on_completed)
        self._worker.cancelled.connect(self._on_cancelled)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _gpu_toggled(self, checked):
        self._gpu_row.setVisible(checked)
        # The panel's height was frozen at whatever sizeHint said when it was expanded, so
        # a row appearing afterwards would be clipped unless we refit.
        if self._options.shown:
            self._options.resize_panel(True)

    def _selected_gpu_id(self):
        """The chosen CUDA index, or None when there is nothing to choose."""
        if not self._gpu_check.isChecked() or not self._gpus:
            return None
        return self._gpu_combo.currentData()

    def _cancel_clicked(self):
        if self._worker is None or not self._worker.isRunning():
            return
        self._worker.cancel()
        self._cancel_button.setEnabled(False)
        self._cancel_button.setText("Cancelling...")
        self._status_label.setText(
            "Cancelling after the current batch is processed.")

    def _set_running(self, running):
        self._run_button.setEnabled(not running)
        self._run_locscale_button.setEnabled(not running)
        self._cancel_button.setEnabled(running)
        self._cancel_button.setText("Cancel")
        self._progress_bar.setVisible(running)
        if running:
            self._progress_bar.setRange(0, 0)      # busy until a stage reports counts

    def _on_status(self, message):
        # Banner rules / blank lines are decoration: log them verbatim (no timestamp) and
        # leave the compact status widget on the last real message.
        stripped = message.strip()
        is_decoration = (stripped == "" or set(stripped) <= set("=-_ "))
        if not is_decoration:
            self._status_label.setText(message)
            from time import strftime
            self.session.logger.info("[{}] {}".format(strftime("%H:%M:%S"), message))
        else:
            self.session.logger.info(message)

    def _on_progress(self, stage, done, total):
        self._progress_bar.setRange(0, total)
        self._progress_bar.setValue(done)
        self._progress_bar.setFormat(f"{stage}: %v/%m")

    def _on_completed(self, results):
        from .cmd import show_locscale_result, show_results
        self._set_running(False)
        self._status_label.setText("Done.")
        # UI thread: safe to touch models
        if self._result_kind == "amplitude_scaling":
            show_locscale_result(self.session, results, self._template)
            reproduced = results.get("locscale")
        else:
            show_results(self.session, results, self._template)
            reproduced = results.get("feature_enhanced")
        self._maybe_report_verification(reproduced)
        self.session.logger.info(
            "LocScale-FEM: if this is useful in your work, please cite "
            "<a href='https://doi.org/{doi}'>doi:{doi}</a>.".format(doi=self._CITATION_DOI),
            is_html=True)

    def _maybe_report_verification(self, reproduced_array):
        """After a verify-seeded run, log the real-space CC vs the published maps."""
        ref = self._verify_reference
        if not ref or reproduced_array is None:
            return
        reproduced = np.asarray(reproduced_array, dtype=np.float32)
        flat = reproduced.ravel()
        best = None
        for vol in ref.get("published", []):
            if getattr(vol, "deleted", False):
                continue
            try:
                other = np.asarray(vol.data.full_matrix(), dtype=np.float32)
            except Exception:
                continue
            if other.shape != reproduced.shape:
                continue
            cc = float(np.corrcoef(flat, other.ravel())[0, 1])
            if best is None or cc > best[1]:
                best = (vol.name, cc)
        if best is None:
            self.session.logger.info(
                "Verification: no published map on the same grid to compare against.")
            return
        bar = "=" * 52
        self.session.logger.info(
            "\n{bar}\n  Verification vs published — EMD-{tok}\n{bar}\n"
            "  Reproduced FEM vs published '{name}'\n"
            "  Real-space cross-correlation (CC) (best) = {cc:.4f}\n{bar}".format(
                bar=bar, tok=ref.get("token"), name=best[0], cc=best[1]))

    def _on_cancelled(self):
        self._set_running(False)
        self._status_label.setText("Cancelled; no maps were opened.")
        self.session.logger.info("LocScale2: cancelled by the user.")

    def _on_failed(self, message):
        self._set_running(False)
        self._status_label.setText("<font color='red'>Failed; see the log.</font>")
        self.session.logger.error("LocScale2 failed.\n" + message)

    def delete(self):
        self._clear_noise_surfaces()
        if self._worker is not None and self._worker.isRunning():
            # Ask first, then wait: otherwise closing the tool blocks ChimeraX for the rest
            # of the run.
            self._worker.cancel()
            self._worker.wait()
        super().delete()
