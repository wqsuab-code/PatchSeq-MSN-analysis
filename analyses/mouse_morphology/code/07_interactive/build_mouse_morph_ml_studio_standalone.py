"""Bundle the Mouse morphology ML studio into one file://-safe HTML document."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "interactive/mouse-morph-ml-studio/dist"
OUTPUT = ROOT / "interactive/mouse-morph-ml-studio/Mouse_M1-M4_ML_validation_studio.html"


def main() -> None:
    html = (DIST / "index.html").read_text(encoding="utf-8")
    css = (DIST / "ml-editor.css").read_text(encoding="utf-8")
    data = (DIST / "ml-data.js").read_text(encoding="utf-8")
    app = (DIST / "ml-editor.js").read_text(encoding="utf-8")
    html = html.replace('<link rel="stylesheet" href="ml-editor.css">', f"<style>{css}</style>")
    html = html.replace('<script src="ml-data.js?v=3"></script><script src="ml-editor.js?v=3"></script>',
                        f"<script>{data}</script><script>{app}</script>")
    OUTPUT.write_text(html, encoding="utf-8")
    print(f"Wrote {OUTPUT} ({OUTPUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
