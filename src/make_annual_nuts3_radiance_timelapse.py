#!/usr/bin/env python3
"""Create an 11-second talk MP4 from locked annual NUTS3 radiance map PNGs.

Each 2014–2024 map is displayed for EXACTLY 1.0 s (30 frames at 30 fps).
The PNGs are outputs of plot_publication_annual_nuts3_radiance_maps.py:
they share geographic extent, colour normalization, year labels and legend.

This is video encoding ONLY; it does not re-evaluate VIIRS, change radiance
values, interpolate images, or replace any PNG/SVG/PDF source file.
"""

from __future__ import annotations

import argparse
import json
import shutil
import struct
import subprocess
from pathlib import Path

from project_config import PROJECT_ROOT as ROOT

YEARS = tuple(range(2014, 2025))
DEFAULT_DIR = ROOT / "figures/publication/annual_nuts3_radiance_2014_2024"
DEFAULT_FILENAME = "annual_nuts3_radiance_timelapse_2014_2024.mp4"
FPS = 30
FRAMES_PER_YEAR = FPS
EXPECTED_FRAMES = len(YEARS) * FRAMES_PER_YEAR
DIMENSIONS = (1920, 1080)


def resolve(path: Path | str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def png_dimensions(path: Path) -> tuple[int, int]:
    """Read PNG IHDR without a third-party image dependency."""
    with path.open("rb") as stream:
        header = stream.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Input is not a PNG: {path}")
    if header[12:16] != b"IHDR":
        raise ValueError(f"PNG header has no IHDR: {path}")
    width, height = struct.unpack(">II", header[16:24])
    if width <= 0 or height <= 0:
        raise ValueError(f"PNG dimensions are invalid: {path}")
    return width, height


def input_inventory(folder: Path) -> list[dict]:
    if not folder.is_dir():
        raise FileNotFoundError(
            f"Annual map output folder missing: {folder}\n"
            "First run: python src/plot_publication_annual_nuts3_radiance_maps.py"
        )
    rows = []
    for year in YEARS:
        path = folder / f"radiance_nuts3_{year}.png"
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing {year} PNG: {path}. Render the annual atlas first."
            )
        width, height = png_dimensions(path)
        rows.append({
            "year": year,
            "filename": path.name,
            "width": width,
            "height": height,
            "bytes": path.stat().st_size,
        })
    if len(rows) != 11 or [r["year"] for r in rows] != list(YEARS):
        raise ValueError("Timelapse must have exactly 11 ordered years 2014–2024")
    return rows


def make_command(ffmpeg: str, folder: Path, output: Path) -> list[str]:
    # -framerate 1 means one sequential source image per second.
    # -r 30 expands each still to 30 identical output frames, with NO
    # morphing/interpolation that might imply an unmeasured annual value.
    # -frames:v fixes an 11-s, 330-frame movie exactly.
    vf = (
        "scale=1920:1080:force_original_aspect_ratio=decrease:flags=lanczos,"
        "pad=1920:1080:(ow-iw)/2:(oh-ih)/2:color=white,"
        "setsar=1"
    )
    return [
        ffmpeg, "-hide_banner", "-loglevel", "warning", "-stats",
        "-y",
        "-framerate", "1",
        "-start_number", "2014",
        "-i", str(folder / "radiance_nuts3_%d.png"),
        "-vf", vf,
        "-r", str(FPS),
        "-frames:v", str(EXPECTED_FRAMES),
        "-c:v", "libx264",
        "-preset", "medium",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        "-movflags", "+faststart",
        str(output),
    ]


def probe_video(output: Path) -> dict:
    ffprobe = shutil.which("ffprobe")
    if not ffprobe:
        return {"ffprobe_available": False, "video_verified": False}
    command = [
        ffprobe, "-v", "error",
        "-select_streams", "v:0",
        "-show_entries", "stream=width,height,nb_frames,r_frame_rate:format=duration",
        "-of", "json",
        str(output),
    ]
    result = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    parsed = json.loads(result.stdout)
    video = parsed.get("streams", [{}])[0]
    fmt = parsed.get("format", {})
    width, height = int(video["width"]), int(video["height"])
    if (width, height) != DIMENSIONS:
        raise ValueError(f"Encoded dimensions unexpected: {(width, height)}")
    rate = str(video.get("r_frame_rate", ""))
    if rate not in {"30/1", "30"}:
        raise ValueError(f"Encoded frame rate unexpected: {rate}")
    frame_count = video.get("nb_frames")
    if frame_count not in (None, "N/A") and int(frame_count) != EXPECTED_FRAMES:
        raise ValueError(
            f"Expected {EXPECTED_FRAMES} video frames, got {frame_count}"
        )
    seconds = float(fmt["duration"])
    if abs(seconds - 11.0) > .06:
        raise ValueError(
            f"Expected 11.0-s MP4 for eleven one-second stills, got {seconds:.3f}s"
        )
    return {
        "ffprobe_available": True,
        "video_verified": True,
        "codec": "H.264",
        "width": width,
        "height": height,
        "frame_rate": rate,
        "frames": int(frame_count) if frame_count not in (None, "N/A") else None,
        "duration_seconds": seconds,
    }


def self_test() -> None:
    assert len(YEARS) == 11 and YEARS[0] == 2014 and YEARS[-1] == 2024
    assert FPS == 30 and EXPECTED_FRAMES == 330
    command = make_command("ffmpeg", Path("/tmp/maps"), Path("/tmp/output.mp4"))
    assert "-framerate" in command and command[command.index("-framerate")+1] == "1"
    assert "-start_number" in command and command[command.index("-start_number")+1] == "2014"
    assert "-frames:v" in command and command[command.index("-frames:v")+1] == "330"
    assert "libx264" in command and "yuv420p" in command
    print("2014–2024 timelapse 11 × 1 second / 1080p output self-test: PASS")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maps-dir", default=DEFAULT_DIR, type=Path,
                    help="Annual radiance PNG folder")
    ap.add_argument("--output", default=None, type=Path,
                    help="Destination MP4 (default alongside source PNGs)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    self_test()
    if args.self_test:
        return

    folder = resolve(args.maps_dir)
    rows = input_inventory(folder)
    output = (
        resolve(args.output)
        if args.output is not None else folder / DEFAULT_FILENAME
    )
    if output.suffix.lower() != ".mp4":
        ap.error("The output filename must end in .mp4")
    output.parent.mkdir(parents=True, exist_ok=True)

    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise FileNotFoundError(
            "ffmpeg is not installed on EcoVis. Install it using "
            "sudo apt update && sudo apt install ffmpeg"
        )

    print("Input: 11 yearly radiance maps, 2014–2024", flush=True)
    for row in rows:
        print(f"  {row['year']}: {row['filename']} "
              f"({row['width']}×{row['height']})", flush=True)
    print("Encoding: 1.0 s/year, 30 fps, 1920×1080 H.264 MP4; "
          "no image interpolation", flush=True)
    cmd = make_command(ffmpeg, folder, output)
    subprocess.run(cmd, check=True)
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"MP4 was not created: {output}")
    video_qa = probe_video(output)

    report = {
        "status": "encoded",
        "input_source": "approved annual 2014–2024 radiance map PNGs",
        "years": list(YEARS),
        "seconds_per_year": 1.0,
        "expected_video_duration_seconds": 11,
        "frames_per_second": FPS,
        "frames_per_year": FRAMES_PER_YEAR,
        "expected_total_frames": EXPECTED_FRAMES,
        "resolution": list(DIMENSIONS),
        "codec": "libx264 H.264",
        "no_tweening_or_interpolation": True,
        "fixed_colours_and_geography_inherited_from_maps": True,
        "existing_PNGs_SVGs_PDFs_modified": False,
        "input_maps": rows,
        "output": str(output),
        "video_probe": video_qa,
    }
    qa_path = output.with_suffix(".json")
    qa_path.write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print("\nANNUAL RADIANCE TIMELAPSE OUTPUT", flush=True)
    print(f"  {output} ({output.stat().st_size/1048576:.2f} MiB)", flush=True)
    print(f"  {qa_path} ({qa_path.stat().st_size/1024:.1f} KiB)", flush=True)
    print("  MP4 QA:", video_qa, flush=True)


if __name__ == "__main__":
    main()
