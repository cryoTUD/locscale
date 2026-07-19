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
                          QPushButton, QSpinBox)

with open(os.path.join(os.path.dirname(__file__), "data", "help_info.json")) as _f:
    help_info = json.load(_f)


class PipelineWorker(QThread):
    """Runs the pipeline off the UI thread.

    Nothing here may touch ChimeraX models: results come back through `finished` and are
    turned into volumes on the UI thread.
    """

    status = Signal(str)
    progress = Signal(str, int, int)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, kwargs):
        super().__init__()
        self._kwargs = kwargs

    def run(self):
        from .pipeline import run_feature_enhance
        try:
            results = run_feature_enhance(
                status_callback=self.status.emit,
                progress_callback=lambda stage, i, n: self.progress.emit(stage, i, n),
                **self._kwargs)
        except Exception as exc:
            import traceback
            self.failed.emit("{}: {}\n\n{}".format(type(exc).__name__, exc, traceback.format_exc()))
            return
        self.finished.emit(results)


class LocScale2Tool(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "LocScale2"
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

        self._map_menu = ModelMenuButton(self.session, class_filter=Volume)
        panel.addWidget(self._row(frame, "Input map:", help_info["input_map_help"], self._map_menu))

        self._mask_menu = ModelMenuButton(self.session, class_filter=Volume)
        panel.addWidget(self._row(frame, "Mask (optional):", help_info["mask_help"], self._mask_menu))

        note = QLabel("<i>If no mask is given, an FDR mask is computed and returned.</i>", frame)
        note.setWordWrap(True)
        panel.addWidget(note)
        return frame

    def _options_panel(self, parent):
        panel = CollapsiblePanel(parent, title="Options")
        frame = panel.content_area
        options = vertical_layout(frame, margins=(0, 0, 0, 0))

        self._model_combo = QComboBox(frame)
        from .emmernet import available_models
        self._model_combo.addItems(available_models())
        options.addWidget(self._row(frame, "Model:", help_info["model_help"], self._model_combo))

        self._mc_spin = QSpinBox(frame); self._mc_spin.setRange(2, 100); self._mc_spin.setValue(15)
        options.addWidget(self._row(frame, "Monte-Carlo iterations:",
                                    help_info["monte_carlo_help"], self._mc_spin))

        self._batch_spin = QSpinBox(frame); self._batch_spin.setRange(1, 128); self._batch_spin.setValue(8)
        options.addWidget(self._row(frame, "Batch size:", help_info["batch_size_help"], self._batch_spin))

        self._window_spin = QSpinBox(frame); self._window_spin.setRange(9, 101); self._window_spin.setValue(25)
        options.addWidget(self._row(frame, "Scaling window:", help_info["window_help"], self._window_spin))

        self._chunk_spin = QSpinBox(frame)
        self._chunk_spin.setRange(128, 65536); self._chunk_spin.setSingleStep(512)
        self._chunk_spin.setValue(4096)
        options.addWidget(self._row(frame, "Scaling chunk:", help_info["chunk_help"], self._chunk_spin))

        self._gpu_check = QCheckBox("Use GPU when available", frame)
        self._gpu_check.setChecked(True)
        self._gpu_check.setToolTip(help_info["gpu_help"])
        options.addWidget(self._gpu_check)
        return panel

    def _run_panel(self, parent):
        frame = QFrame(parent)
        panel = vertical_layout(frame, margins=(0, 8, 0, 0))

        self._run_button = QPushButton("Run feature enhancement", frame)
        self._run_button.clicked.connect(self._run_clicked)
        panel.addWidget(self._run_button)

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
        mask = mask_volume.data.full_matrix() if mask_volume is not None else None
        if mask is not None and mask.shape != emmap.shape:
            self._status_label.setText(
                "<font color='red'>Mask shape {} does not match map {}.</font>".format(
                    mask.shape, emmap.shape))
            return

        self._template = input_map
        self._set_running(True)

        self._worker = PipelineWorker(dict(
            emmap=emmap,
            apix=float(input_map.data.step[0]),
            mask=mask,
            model_type=self._model_combo.currentText(),
            monte_carlo_iterations=self._mc_spin.value(),
            batch_size=self._batch_spin.value(),
            window_size=self._window_spin.value(),
            scaling_chunk=self._chunk_spin.value(),
            use_gpu=self._gpu_check.isChecked(),
        ))
        self._worker.status.connect(self._on_status)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _set_running(self, running):
        self._run_button.setEnabled(not running)
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

    def _on_finished(self, results):
        from .cmd import show_results
        self._set_running(False)
        self._status_label.setText("Done.")
        show_results(self.session, results, self._template)   # UI thread: safe to touch models

    def _on_failed(self, message):
        self._set_running(False)
        self._status_label.setText("<font color='red'>Failed &mdash; see the log.</font>")
        self.session.logger.error("LocScale2 failed.\n" + message)

    def delete(self):
        if self._worker is not None and self._worker.isRunning():
            self._worker.wait()
        super().delete()
