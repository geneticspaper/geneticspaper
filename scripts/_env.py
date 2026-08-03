"""Minimal .env loader (standard library only).

Reads key=value pairs from a project-root `.env` file into os.environ without
overriding variables already set in the environment. Keeps local, machine-
specific paths (e.g. the location of the source supplementary workbooks) out of
version control. See `.env.example` for the recognized keys.
"""
from __future__ import annotations

import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENV_PATH = os.path.join(ROOT, ".env")


def load_env(path: str = ENV_PATH) -> None:
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = value


def get(key: str, default: str = "") -> str:
    load_env()
    return os.environ.get(key, default)
