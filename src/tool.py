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

from chimerax.core.tools import ToolInstance
from chimerax.map import Volume
from chimerax.ui import MainToolWindow
from chimerax.ui.widgets import CollapsiblePanel, ModelMenuButton, vertical_layout
from Qt.QtCore import QThread, Signal
from Qt.QtWidgets import (QCheckBox, QComboBox, QFrame, QHBoxLayout, QLabel, QProgressBar,
                          QPushButton, QSpinBox, QLineEdit, QVBoxLayout)


from .utils import (
    round_up_proper,
    round_up_to_even,
)

with open(os.path.join(os.path.dirname(__file__), "data", "help_info.json")) as _f:
    help_info = json.load(_f)


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

    def __init__(self, kwargs):
        super().__init__()
        self._kwargs = kwargs
        self._cancel_requested = False

    def cancel(self):
        """Ask the run to stop. Safe to call from the UI thread: it only sets a flag."""
        self._cancel_requested = True

    def run(self):
        from .pipeline import Cancelled, run_feature_enhance

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
            results = run_feature_enhance(
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

        self.tool_window = MainToolWindow(self)
        parent = self.tool_window.ui_area
        layout = vertical_layout(parent, margins=(5, 5, 5, 5))
        layout.addWidget(self._inputs_panel(parent))
        layout.addWidget(self._options_panel(parent))
        layout.addWidget(self._run_panel(parent))
        layout.addStretch(1)
        self.tool_window.manage("side")

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

    def _inputs_panel(self, parent):
        frame = QFrame(parent)
        panel = vertical_layout(frame, margins=(0, 0, 0, 8))
        panel.addWidget(QLabel("<b>Inputs</b>", frame))

        # Add input unsharpened map
        self._map_menu = ModelMenuButton(self.session, class_filter=Volume)
        panel.addWidget(self._row(frame, "Input map:", help_info["input_map_help"], self._map_menu))
        # Add optional mask
        self._mask_menu = ModelMenuButton(
            self.session, 
            class_filter=Volume, 
            no_value_button_text="No model chosen",
            no_value_menu_text="None",
            autoselect="none",
        )
        panel.addWidget(self._row(frame, "Mask (optional):", help_info["mask_help"], self._mask_menu))
        note = QLabel("<i>If no mask is given, an FDR mask is computed and returned.</i>", frame)
        note.setWordWrap(True)
        panel.addWidget(note)

        # Add point group symmetry as text input 
        self._point_group_symmetry_menu = QLineEdit(frame)
        panel.addWidget(self._row(frame, "Point group symmetry:", help_info["point_group_symmetry_help"], self._point_group_symmetry_menu))

        # Add helical symmetry checkbox default to unchecked
        self._helical_symmetry_checkbox = QCheckBox("Helical symmetry", frame)
        self._helical_symmetry_checkbox.setToolTip(help_info["helical_symmetry_help"])
        self._helical_symmetry_checkbox.toggled.connect(self._helical_symmetry_toggled)
        self._helical_symmetry_checkbox.setChecked(False)
        panel.addWidget(self._row(frame, "Helical symmetry:", help_info["helical_symmetry_help"], self._helical_symmetry_checkbox))
        
        # Add helical rise/twist inputs. They live in their own panel that is only shown
        # while the helical checkbox is ticked (see _helical_symmetry_toggled).
        self._helical_symmetry_panel = QFrame(frame)
        helical_layout = QVBoxLayout(self._helical_symmetry_panel)
        helical_layout.setContentsMargins(0, 0, 0, 0)
        # input rise in Angstroms
        self._helical_rise_menu = QLineEdit(self._helical_symmetry_panel, placeholderText="0.0")
        helical_layout.addWidget(self._row(self._helical_symmetry_panel, "Rise (A):",
                                           help_info["helical_symmetry_help"], self._helical_rise_menu))
        # input twist in degrees
        self._helical_twist_menu = QLineEdit(self._helical_symmetry_panel, placeholderText="0.0")
        helical_layout.addWidget(self._row(self._helical_symmetry_panel, "Twist (deg):",
                                           help_info["helical_symmetry_help"], self._helical_twist_menu))
        self._helical_symmetry_panel.setVisible(False)
        panel.addWidget(self._helical_symmetry_panel)
        return frame

    def _helical_symmetry_toggled(self, checked):
        """Show the rise/twist inputs only while helical symmetry is requested."""
        self._helical_symmetry_panel.setVisible(checked)

    def _options_panel(self, parent):
        self._options = panel = CollapsiblePanel(parent, title="Advanced Options")
        frame = panel.content_area
        # CollapsiblePanel already installs a QVBoxLayout on content_area. Calling
        # vertical_layout() here would try to set a second layout on the same widget, which
        # Qt refuses: the rows below would never be laid out, content_area.sizeHint() would
        # stay at zero, and expanding the panel would reveal nothing.
        options = frame.layout()

        self._model_combo = QComboBox(frame)
        from .emmernet import available_models
        self._model_combo.addItems(available_models())
        options.addWidget(self._row(frame, "Model:", help_info["model_help"], self._model_combo))

        
        self._mc_spin = QSpinBox(frame); self._mc_spin.setRange(1, 100); self._mc_spin.setValue(15)
        options.addWidget(self._row(frame, "Monte-Carlo iterations:",
                                    help_info["monte_carlo_help"], self._mc_spin))

        self._batch_spin = QSpinBox(frame); self._batch_spin.setRange(1, 64); self._batch_spin.setValue(8)
        options.addWidget(self._row(frame, "Batch size:", help_info["batch_size_help"], self._batch_spin))

        # Cube size is deliberately not exposed: EMmerNet's padding and crops assume 32.
        self._stride_spin = QSpinBox(frame); self._stride_spin.setRange(1, 32); self._stride_spin.setValue(16)
        options.addWidget(self._row(frame, "Cube stride:", help_info["stride_help"], self._stride_spin))

        self._window_spin = QSpinBox(frame); self._window_spin.setRange(11, 45); self._window_spin.setValue(25)
        options.addWidget(self._row(frame, "Scaling window:", help_info["window_help"], self._window_spin))

        self._chunk_spin = QSpinBox(frame)
        self._chunk_spin.setRange(128, 65536); self._chunk_spin.setSingleStep(512)
        self._chunk_spin.setValue(4096)
        options.addWidget(self._row(frame, "Scaling chunk:", help_info["chunk_help"], self._chunk_spin))

        self._gpu_check = QCheckBox("Use GPU when available", frame)
        self._gpu_check.setChecked(True)
        self._gpu_check.setToolTip(help_info["gpu_help"])
        self._gpu_check.toggled.connect(self._gpu_toggled)
        options.addWidget(self._gpu_check)

        from .emmernet import available_gpus
        self._gpus = available_gpus()
        self._gpu_combo = QComboBox(frame)
        for index, name in self._gpus:
            self._gpu_combo.addItem("{}: {}".format(index, name), index)
        if not self._gpus:
            # MPS and CPU have no device index to pick, so leave the control in place but
            # inert rather than implying a choice that does not exist.
            self._gpu_combo.addItem("Using CPU", None)
            self._gpu_combo.setEnabled(False)
        self._gpu_row = self._row(frame, "GPU:", help_info["gpu_id_help"], self._gpu_combo)
        options.addWidget(self._gpu_row)
        return panel

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
        input_map = self._map_menu.value
        mask_volume = self._mask_menu.value

        if input_map is None:
            self._status_label.setText("<font color='red'>Select an input map.</font>")
            return
        if mask_volume is input_map:
            self._status_label.setText("<font color='red'>Map and mask are the same volume.</font>")
            return

        emmap = input_map.data.full_matrix()
        apix = float(input_map.data.step[0])

        window_size_angstroms = self._window_spin.value()
        window_size_pix = round_up_to_even(window_size_angstroms / apix)

        mask = mask_volume.data.full_matrix() if mask_volume is not None else None
        if mask is not None and mask.shape != emmap.shape:
            self._status_label.setText(
                "<font color='red'>Mask shape {} does not match map {}.</font>".format(
                    mask.shape, emmap.shape))
            return

        self._template = input_map
        self._set_running(True)

        # Get the point group symmetry and helical symmetry parameters.
        # An empty field means "no symmetry": placeholderText is only a hint, not a value,
        # so .text() is "" when the user leaves it blank -- default that to C1.
        pg = self._point_group_symmetry_menu.text().strip() or "C1"
        helical_symmetry = self._helical_symmetry_checkbox.isChecked()
        rise = float(self._helical_rise_menu.text()) if helical_symmetry else None
        twist = float(self._helical_twist_menu.text()) if helical_symmetry else None

        self._worker = PipelineWorker(dict(
            emmap=emmap,
            apix=apix,
            mask=mask,
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
        self._cancel_button.setEnabled(running)
        self._cancel_button.setText("Cancel")
        self._progress_bar.setVisible(running)
        if running:
            self._progress_bar.setRange(0, 0)      # busy until a stage reports counts

    def _on_status(self, message):
        self._status_label.setText(message)
        self.session.logger.info("LocScale2: " + message)

    def _on_progress(self, stage, done, total):
        self._progress_bar.setRange(0, total)
        self._progress_bar.setValue(done)
        self._progress_bar.setFormat(f"{stage}: %v/%m")

    def _on_completed(self, results):
        from .cmd import show_results
        self._set_running(False)
        self._status_label.setText("Done.")
        show_results(self.session, results, self._template)   # UI thread: safe to touch models

    def _on_cancelled(self):
        self._set_running(False)
        self._status_label.setText("Cancelled; no maps were opened.")
        self.session.logger.info("LocScale2: cancelled by the user.")

    def _on_failed(self, message):
        self._set_running(False)
        self._status_label.setText("<font color='red'>Failed; see the log.</font>")
        self.session.logger.error("LocScale2 failed.\n" + message)

    def delete(self):
        if self._worker is not None and self._worker.isRunning():
            # Ask first, then wait: otherwise closing the tool blocks ChimeraX for the rest
            # of the run.
            self._worker.cancel()
            self._worker.wait()
        super().delete()
