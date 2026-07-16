"""Host-to-container path translation for a bounded Docker workspace."""

from pathlib import Path, PurePosixPath


class DockerWorkspace:
    """Map paths under one host directory into one container directory."""

    def __init__(self, host_root, container_root):
        self.host_root = Path(host_root).expanduser().resolve()
        self.container_root = PurePosixPath(str(container_root))
        if not self.container_root.is_absolute():
            raise ValueError("Docker container workspace must be absolute")
        if self.container_root == PurePosixPath('/'):
            raise ValueError("Docker container workspace cannot be root")
        if '..' in self.container_root.parts:
            raise ValueError("Docker container workspace cannot contain '..'")

    @property
    def volume_spec(self):
        """Return a Docker ``-v`` bind specification."""
        return f"{self.host_root}:{self.container_root}"

    def container_path(self, host_path):
        """Translate a workspace path to its container path.

        Relative paths are interpreted from the current directory when that
        resolves inside the workspace, and from ``host_root`` otherwise.
        """
        candidate = Path(host_path).expanduser()
        if candidate.is_absolute():
            resolved = candidate.resolve()
        else:
            cwd_candidate = (Path.cwd() / candidate).resolve()
            if self._is_within(cwd_candidate, self.host_root):
                resolved = cwd_candidate
            else:
                resolved = (self.host_root / candidate).resolve()

        try:
            relative = resolved.relative_to(self.host_root)
        except ValueError as exc:
            raise ValueError(
                f"Docker path is outside the mounted workspace: {resolved}; "
                f"workspace={self.host_root}"
            ) from exc

        return str(self.container_root.joinpath(*relative.parts))

    @staticmethod
    def _is_within(path, root):
        try:
            path.relative_to(root)
        except ValueError:
            return False
        return True
