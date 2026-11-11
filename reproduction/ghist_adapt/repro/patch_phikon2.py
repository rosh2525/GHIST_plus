"""Patch a GHIST+ checkout so the image encoder can be swapped to the open owkin/phikon-v2 (GHIST_ENCODER=phikon2).
Default behaviour (UNI2-h) is unchanged. Usage: python patch_phikon2.py /path/to/GHIST_plus
"""
import sys
from pathlib import Path

ENC = '''

class PhikonV2Encoder(nn.Module):
    """Open replacement for UNI2-h (used when GHIST_ENCODER=phikon2): owkin/phikon-v2, a DINOv2 ViT-L/16
    pathology encoder (1024-d), loaded with transformers. Same output contract as Uni2HEncoder:
    (B, C, H/16, W/16) patch-token features. ImageNet input normalisation, as for UNI2-h."""
    patch_size = 16

    def __init__(self, pretrained: bool = True):
        super().__init__()
        from transformers import AutoConfig, AutoModel
        if pretrained:
            self.vit = AutoModel.from_pretrained("owkin/phikon-v2")
        else:
            self.vit = AutoModel.from_config(AutoConfig.from_pretrained("owkin/phikon-v2"))
        self.out_channels = self.vit.config.hidden_size      # 1024

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, _, H, W = x.shape
        seq = self.vit(pixel_values=x).last_hidden_state    # (B, 1 + N, C)
        h, w = H // self.patch_size, W // self.patch_size
        tokens = seq[:, -h * w:, :]
        return tokens.transpose(1, 2).reshape(B, -1, h, w)

'''

p = Path(sys.argv[1]) / "model/backbone.py"
s = p.read_text()
if "PhikonV2Encoder" in s:
    sys.exit("already patched")
s = s.replace("class Uni2HAdapter(nn.Module):", ENC + "class Uni2HAdapter(nn.Module):", 1)
old = "        self.enc = Uni2HEncoder(pretrained=pretrained)   # ViT-H/14"
assert old in s
s = s.replace(old, '        import os\n        self.enc = (PhikonV2Encoder(pretrained=pretrained) if os.getenv("GHIST_ENCODER", "uni2h") == "phikon2"\n'
                   '                    else Uni2HEncoder(pretrained=pretrained))   # ViT-H/14 by default', 1)
old2 = "        pad_h = (14 - H % 14) % 14\n        pad_w = (14 - W % 14) % 14"
assert old2 in s
s = s.replace(old2, '        ps = getattr(self.enc, "patch_size", 14)\n        pad_h = (ps - H % ps) % ps\n        pad_w = (ps - W % ps) % ps', 1)
p.write_text(s)
print("backbone patched:", p)
