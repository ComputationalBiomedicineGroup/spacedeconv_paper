#!/usr/bin/env bash
set -euo pipefail

figures_dir="./export/figures"
top="./data/spacedeconv_overview.png"
bottom="${figures_dir}/fig_1_brain_2row.png"
output="${figures_dir}/fig_1.png"
gap_px=550
label_font="DejaVu-Sans-Bold"
label_pointsize=200

tmp_top="$(mktemp "${figures_dir}/fig_1_overview_resized.XXXXXX.png")"
tmp_gap="$(mktemp "${figures_dir}/fig_gap.XXXXXX.png")"
trap 'rm -f "${tmp_top}" "${tmp_gap}"' EXIT

target_width="$(identify -format '%w' "${bottom}")"

convert "${top}" \
  -resize "${target_width}x" \
  -background white \
  -alpha remove \
  -alpha off \
  "${tmp_top}"

convert "${tmp_top}" \
  -font "${label_font}" \
  -pointsize "${label_pointsize}" \
  -fill black \
  -gravity northwest \
  -annotate +0+0 "A" \
  "${tmp_top}"

convert -size "${target_width}x${gap_px}" xc:white "${tmp_gap}"

convert "${tmp_top}" "${tmp_gap}" "${bottom}" \
  -background white \
  -gravity center \
  -append \
  "${output}"

identify "${output}"
