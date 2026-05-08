"""Tests for the ENSG → gene symbol lookup (T-13 support)."""

from __future__ import annotations

from pathlib import Path

from genomescope.data import symbols


def test_fetch_symbols_uses_cache(tmp_path: Path, monkeypatch) -> None:
    cache_file = tmp_path / "ensg_to_symbol.json"
    monkeypatch.setattr(symbols, "CACHE", cache_file)

    # pre-seed cache so no network calls
    cache_file.write_text('{"ENSG00000001": "APP", "ENSG00000002": "SOD1"}')

    # monkey-patch requests so any attempted network call fails loudly
    def boom(*args, **kwargs):
        raise AssertionError("should not hit network when cache is warm")
    monkeypatch.setattr(symbols.requests, "post", boom)

    out = symbols.fetch_symbols(["ENSG00000001.3", "ENSG00000002"])
    assert out == {"ENSG00000001": "APP", "ENSG00000002": "SOD1"}


def test_fetch_symbols_handles_network_failure(tmp_path: Path, monkeypatch) -> None:
    """Offline or rate-limited: fall back to self-mapping."""
    cache_file = tmp_path / "ensg_to_symbol.json"
    monkeypatch.setattr(symbols, "CACHE", cache_file)

    class BoomResponse:
        def raise_for_status(self) -> None:
            raise symbols.requests.ConnectionError("offline")

    def offline(*args, **kwargs):
        raise symbols.requests.ConnectionError("offline")
    monkeypatch.setattr(symbols.requests, "post", offline)

    out = symbols.fetch_symbols(["ENSG00000003"])
    # couldn't resolve → mapped to itself
    assert out == {"ENSG00000003": "ENSG00000003"}


def test_fetch_symbols_strips_version(tmp_path: Path, monkeypatch) -> None:
    cache_file = tmp_path / "ensg_to_symbol.json"
    monkeypatch.setattr(symbols, "CACHE", cache_file)
    cache_file.write_text('{"ENSG00000001": "APP"}')
    out = symbols.fetch_symbols(["ENSG00000001.7"])
    assert out == {"ENSG00000001": "APP"}
