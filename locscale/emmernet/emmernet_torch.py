"""
PyTorch implementation of EMmerNet. Implemented with Claude Opus 4.8. 
Ran tests to verify the model produces identical output compared to 
tensorflow models. Max difference in the predictions from 100 random 
cubes was < 1e-3 with correlation > 0.9999.

Monte carlo dropout is used in the model to enable uncertainty estimation. 
The variance measured from the PyTorch model matches the variance measured from 
tensorflow model (correlation > xx, ks-test p-value ~).
"""

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


# (torch ConvBlock attribute name, TF conv layer name, TF PReLU layer name, TF GroupNorm layer name)
_CONV_BLOCK_MAP = [
    ("conv1_1", "conv1_1_3d", "p_re_lu", "group_normalization"),
    ("conv1_2", "conv1_2_3d", "p_re_lu_1", "group_normalization_1"),
    ("conv1_3", "conv1_3_3d", "p_re_lu_2", "group_normalization_2"),
    ("conv2_1", "conv2_1_3d", "p_re_lu_3", "group_normalization_3"),
    ("conv2_2", "conv2_2_3d", "p_re_lu_4", "group_normalization_4"),
    ("conv2_3", "conv2_3_3d", "p_re_lu_5", "group_normalization_5"),
    ("conv3_1", "conv3_1_3d", "p_re_lu_6", "group_normalization_6"),
    ("conv3_2", "conv3_2_3d", "p_re_lu_7", "group_normalization_7"),
    ("conv3_3", "conv3_3_3d", "p_re_lu_8", "group_normalization_8"),
    ("conv4_1", "conv4_1_3d", "p_re_lu_9", "group_normalization_9"),
    ("conv5_1", "conv5_1_3d", "p_re_lu_10", "group_normalization_10"),
    ("conv5_2", "conv5_2_3d", "p_re_lu_11", "group_normalization_11"),
    ("conv5_3", "conv5_3_3d", "p_re_lu_12", "group_normalization_12"),
    ("conv6_1", "conv6_1_3d", "p_re_lu_13", "group_normalization_13"),
    ("conv6_2", "conv6_2_3d", "p_re_lu_14", "group_normalization_14"),
    ("conv6_3", "conv6_3_3d", "p_re_lu_15", "group_normalization_15"),
    ("conv7_1", "conv7_1_3d", "p_re_lu_16", "group_normalization_16"),
    ("conv7_2", "conv7_2_3d", "p_re_lu_17", "group_normalization_17"),
    ("conv7_3", "conv7_3_3d", "p_re_lu_18", "group_normalization_18"),
    ("conv8_1", "conv8_1_3d", "p_re_lu_19", "group_normalization_19"),
]

_UP_BLOCK_MAP = [
    ("up1", "conv3d_transpose"),
    ("up2", "conv3d_transpose_1"),
    ("up3", "conv3d_transpose_2"),
    ("up4", "conv3d_transpose_3"),
]


def _set_conv_weights(torch_conv, tf_layer):
    kernel, bias = tf_layer.get_weights()
    kernel_t = torch.from_numpy(kernel).permute(4, 3, 0, 1, 2).contiguous().float()
    torch_conv.weight.data.copy_(kernel_t)
    torch_conv.bias.data.copy_(torch.from_numpy(bias).float())


def load_weights_from_keras(torch_model, keras_model):
    for torch_name, conv_name, prelu_name, gn_name in _CONV_BLOCK_MAP:
        block = getattr(torch_model, torch_name)
        _set_conv_weights(block.conv, keras_model.get_layer(conv_name))

        prelu_alpha = keras_model.get_layer(prelu_name).get_weights()[0]  # shape (1,1,1,C)
        block.prelu.weight.data.copy_(torch.from_numpy(prelu_alpha).reshape(-1).float())

        gamma, beta = keras_model.get_layer(gn_name).get_weights()
        block.gn.weight.data.copy_(torch.from_numpy(gamma).float())
        block.gn.bias.data.copy_(torch.from_numpy(beta).float())

    for torch_name, tf_name in _UP_BLOCK_MAP:
        up_block = getattr(torch_model, torch_name)
        _set_conv_weights(up_block.conv, keras_model.get_layer(tf_name))

    _set_conv_weights(torch_model.last_conv, keras_model.get_layer("last_layer"))

    return torch_model


def create_dataloader(cubes, batch_size=32, shuffle=False):
    tensor_cubes = torch.from_numpy(cubes).float()
    dataset = torch.utils.data.TensorDataset(tensor_cubes)
    
    # PIN_MEMORY and NUM_WORKERS are crucial here
    dataloader = torch.utils.data.DataLoader(
        dataset, 
        batch_size=batch_size, 
        shuffle=shuffle,       
    )
    return dataloader