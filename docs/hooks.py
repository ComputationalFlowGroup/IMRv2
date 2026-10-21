"""Use the repository README as the local MkDocs home without editing it."""

from pathlib import Path
import re
from urllib.parse import urlsplit


def on_page_markdown(markdown, page, config, files):
    if page.file.src_uri != "index.md":
        return markdown

    readme = Path(config.config_file_path).parent / "README.md"
    markdown = readme.read_text(encoding="utf-8")
    repository = config["repo_url"].rstrip("/")

    def destination(match):
        target = match.group("target")
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or target.startswith(("#", "/")):
            return match.group(0)
        if target.startswith("docs/"):
            # MkDocs resolves this relative to docs/index.md and validates it.
            target = target.removeprefix("docs/")
        else:
            # Files outside docs must resolve on GitHub, including LICENSE.
            target = f"{repository}/blob/main/{target}"
        return f"{match.group('prefix')}{target}{match.group('suffix')}"

    return re.sub(
        r"(?P<prefix>\]\()(?P<target>[^\s)]+)(?P<suffix>\))",
        destination,
        markdown,
    )
