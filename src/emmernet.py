"""EMmerNet: the 3D U-Net, and Monte-Carlo dropout inference.

The architecture is the hand-written PyTorch port of the original TensorFlow EMmerNet,
validated against it to correlation > 0.9999 (max abs difference < 1e-3 over 100 random
cubes) with the Monte-Carlo dropout distributions statistically indistinguishable.

Only the high-context checkpoint ships with the bundle; see available_models().
"""
import os

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class ConvBlock(nn.Module):
    def __init__(self, in_ch, out_ch, stride=1):
        super().__init__()
        if stride == 1:
            # TF 'same' padding for k=5, s=1 is exactly symmetric: pad=2 each side.
            self.pre_pad = None
            self.conv = nn.Conv3d(in_ch, out_ch, kernel_size=5, stride=1, padding=2)
        else:
            # TF 'same' padding for k=5, s=2 on our (32/16/8) inputs is asymmetric: 1 before, 2 after.
            self.pre_pad = (1, 2, 1, 2, 1, 2)  # F.pad order: (W_l,W_r, H_l,H_r, D_l,D_r)
            self.conv = nn.Conv3d(in_ch, out_ch, kernel_size=5, stride=2, padding=0)
        self.prelu = nn.PReLU(num_parameters=out_ch)
        self.gn = nn.GroupNorm(num_groups=8, num_channels=out_ch, eps=0.001)
        self.dropout = nn.Dropout(p=0.5)

    def forward(self, x):
        if self.pre_pad is not None:
            x = F.pad(x, self.pre_pad)
        x = self.conv(x)
        x = self.prelu(x)
        x = self.gn(x)
        x = self.dropout(x)
        return x

class UpBlock(nn.Module):
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.conv = nn.ConvTranspose3d(in_ch, out_ch, kernel_size=5, stride=2, padding=0, output_padding=0)

    def forward(self, x):
        x = self.conv(x)
        return x[:, :, 1:-2, 1:-2, 1:-2]

class EMmerNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1_1 = ConvBlock(1, 32, stride=1)
        self.conv1_2 = ConvBlock(32, 32, stride=1)
        self.conv1_3 = ConvBlock(32, 32, stride=2)

        self.conv2_1 = ConvBlock(32, 64, stride=1)
        self.conv2_2 = ConvBlock(64, 64, stride=1)
        self.conv2_3 = ConvBlock(64, 64, stride=2)

        self.conv3_1 = ConvBlock(64, 128, stride=1)
        self.conv3_2 = ConvBlock(128, 128, stride=1)
        self.conv3_3 = ConvBlock(128, 128, stride=2)

        self.conv4_1 = ConvBlock(128, 128, stride=1)
        self.up1 = UpBlock(128, 128)

        self.conv5_1 = ConvBlock(256, 128, stride=1)
        self.conv5_2 = ConvBlock(128, 128, stride=1)
        self.conv5_3 = ConvBlock(128, 128, stride=1)
        self.up2 = UpBlock(128, 64)

        self.conv6_1 = ConvBlock(128, 64, stride=1)
        self.conv6_2 = ConvBlock(64, 64, stride=1)
        self.conv6_3 = ConvBlock(64, 64, stride=1)
        self.up3 = UpBlock(64, 32)

        self.conv7_1 = ConvBlock(64, 32, stride=1)
        self.conv7_2 = ConvBlock(32, 32, stride=1)
        self.conv7_3 = ConvBlock(32, 32, stride=1)
        self.up4 = UpBlock(32, 16)

        self.conv8_1 = ConvBlock(16, 8, stride=2)
        self.last_conv = nn.Conv3d(8, 1, kernel_size=5, stride=1, padding=2)

    def forward(self, x):
        # x: (N, 1, 32, 32, 32) channels-first
        # Skip connections tap off the FIRST conv block in each encoder stage
        # (verified against the model's actual Concatenate inbound_nodes), not the second --
        # the second block sits on the main path but its own output is never reused.
        e1 = self.conv1_1(x)
        d1 = self.conv1_3(self.conv1_2(e1))
        e2 = self.conv2_1(d1)
        d2 = self.conv2_3(self.conv2_2(e2))
        e3 = self.conv3_1(d2)
        d3 = self.conv3_3(self.conv3_2(e3))

        b = self.conv4_1(d3)
        u1 = self.up1(b)
        c1 = torch.cat([u1, e3], dim=1)
        x5 = self.conv5_3(self.conv5_2(self.conv5_1(c1)))
        u2 = self.up2(x5)
        c2 = torch.cat([u2, e2], dim=1)
        x6 = self.conv6_3(self.conv6_2(self.conv6_1(c2)))
        u3 = self.up3(x6)
        c3 = torch.cat([u3, e1], dim=1)
        x7 = self.conv7_3(self.conv7_2(self.conv7_1(c3)))

        u4 = self.up4(x7)  # no skip connection here -- matches the original's tail quirk
        x8 = self.conv8_1(u4)
        out = self.last_conv(x8)
        return out

CHECKPOINTS = {
    "high_context": "emmernet_highcontext.pt",
    "low_context": "emmernet_lowcontext.pt",
}


def data_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def available_models():
    """Model variants whose weights are actually present.

    Only high_context ships by default: each checkpoint is ~91 MB, and carrying both makes
    the bundle slow to build and awkward to distribute. Drop emmernet_lowcontext.pt into
    data/ and low_context becomes selectable again with no code change.
    """
    present = [name for name, fn in CHECKPOINTS.items()
               if os.path.exists(os.path.join(data_dir(), fn))]
    return present or ["high_context"]      # never return empty; load_emmernet reports the miss


def available_gpus():
    """Selectable GPUs as (index, name) pairs; empty when there is nothing to choose from.

    Only CUDA is enumerated. Apple MPS exposes a single unindexed device, so there is no
    choice to offer there, and CPU obviously has none either.
    """
    if torch.cuda.is_available():
        return [(i, torch.cuda.get_device_name(i)) for i in range(torch.cuda.device_count())]

    elif torch.backends.mps.is_available():
        return [(0, "Apple MPS")]
    
    else:
        return []
    


def get_device(prefer_gpu=True, gpu_id=None):
    """CUDA if available, else Apple MPS, else CPU.

    `gpu_id` picks a specific CUDA device. It is ignored on MPS and CPU, which have no
    device index; an out-of-range index is an error rather than a silent fallback, since
    quietly running on the wrong GPU is worse than refusing to start.
    """
    if not prefer_gpu:
        return torch.device("cpu")
    if torch.cuda.is_available():
        if gpu_id is None:
            return torch.device("cuda")
        count = torch.cuda.device_count()
        if not 0 <= int(gpu_id) < count:
            raise ValueError(
                f"GPU {gpu_id} does not exist; {count} CUDA device(s) present (0-{count - 1})")
        return torch.device("cuda", int(gpu_id))
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_emmernet(model_type="high_context", device=None, data_dir_override=None):
    if model_type not in CHECKPOINTS:
        raise ValueError(f"unknown model_type {model_type!r}; choose from {sorted(CHECKPOINTS)}")
    directory = data_dir_override or data_dir()
    path = os.path.join(directory, CHECKPOINTS[model_type])
    if not os.path.exists(path):
        raise FileNotFoundError(f"EMmerNet weights missing: {path}")

    device = device or get_device()
    model = EMmerNet()
    model.load_state_dict(torch.load(path, map_location="cpu"))
    return model.to(device), device


def _to_nchw(cubes):
    """Return (N, 1, D, H, W) from either cube layout.

    get_cubes() has returned both over time: channels-last (N, D, H, W, 1) in the
    TensorFlow-era code and channels-first (N, 1, D, H, W) in the PyTorch port. Accept
    either rather than depend on which one is current.
    """
    cubes = np.asarray(cubes, dtype=np.float32)
    if cubes.ndim == 4:                       # (N, D, H, W) -- no channel axis
        return cubes[:, None]
    if cubes.ndim != 5:
        raise ValueError(f"expected 4D or 5D cube array, got shape {cubes.shape}")
    if cubes.shape[1] == 1:                   # already (N, 1, D, H, W)
        return cubes
    if cubes.shape[-1] == 1:                  # (N, D, H, W, 1) channels-last
        return np.ascontiguousarray(cubes.transpose(0, 4, 1, 2, 3))
    raise ValueError(f"cannot find the channel axis in cube array of shape {cubes.shape}")


def predict_monte_carlo(cubes, model, device, num_samples=15, batch_size=8,
                        progress_callback=None):
    """Monte-Carlo dropout inference.

    Returns (mean, variance), each (N, cube, cube, cube).

    model.train() enables dropout. EMmerNet has no BatchNorm, and GroupNorm/PReLU are
    unaffected by train/eval mode, so this switches on MC dropout and nothing else.
    """
    model.train()
    cubes = _to_nchw(cubes)
    n = len(cubes)

    means, variances = [], []
    with torch.inference_mode():
        for start in range(0, n, batch_size):
            x = torch.from_numpy(cubes[start:start + batch_size]).to(device)

            samples = torch.stack([model(x) for _ in range(num_samples)], dim=0)
            # unbiased=False matches np.var's default, which the reference used
            var, mean = torch.var_mean(samples, dim=0, unbiased=False)

            means.append(mean.squeeze(1).cpu().numpy())        # (B, 1, D, H, W) -> (B, D, H, W)
            variances.append(var.squeeze(1).cpu().numpy())
            if progress_callback is not None:
                progress_callback(min(start + batch_size, n), n)

    return np.concatenate(means, axis=0), np.concatenate(variances, axis=0)
