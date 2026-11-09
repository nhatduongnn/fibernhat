"""Every skill folder must be registered in CAPABILITIES.md (and have a SKILL.md)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_every_skill_registered():
    caps = (ROOT / "CAPABILITIES.md").read_text()
    for d in (ROOT / ".claude" / "skills").iterdir():
        if d.is_dir() and not d.name.startswith("_"):
            assert (d / "SKILL.md").exists(), d
            assert f"`{d.name}`" in caps, f"{d.name} missing from CAPABILITIES.md"


def test_no_secrets_in_repo():
    import re
    pat = re.compile(r"gh[ops]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}")
    for p in ROOT.rglob("*"):
        if p.is_file() and ".venv" not in p.parts and ".git" not in p.parts and p.suffix in (
                ".py", ".md", ".yaml", ".yml", ".sbatch", ".sh", ".toml", ".cfg", ".txt", ".tsv"):
            assert not pat.search(p.read_text(errors="ignore")), p
