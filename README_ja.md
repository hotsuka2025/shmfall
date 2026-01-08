<table>
	<thead>
		<tr>
			<th style="text-align:center"><a href="README.md">English</a></th>
			<th style="text-align:center">日本語</th>
		</tr>
	</thead>
</table>

# shmfall

**shmfall** は WIN 形式（卜部 1994）の時系列データを、`shmdump` が生成するプレーンテキストストリームから読み取り、リアルタイムでウォーターフォール表示する Python ツールです。現バージョンでは共有メモリを直接扱わず、`shmdump -tq` 等の標準出力をパイプで受け取って動作します。

本 README は shmfall.py の簡単な使用案内です。詳細は下記文献を参照してください。

- 大塚宏徳・田中伸一・篠原雅尚, 2025, WINシステムに対応したリアルタイムWaterfallプロット描画ツール, 東京大学地震研究所技術研究報告, 31号

WINシステム(卜部・束田 1992)については下記を参照してください。

- ["WINシステムとは?"](https://wwweic.eri.u-tokyo.ac.jp/WIN/)

## 概要
- 入力: `shmdump` のプレーンテキストストリーム（または標準入力経由の同等データ）
- 出力: リアルタイムのウォーターフォール表示（Matplotlib を使用）
- 主要実装: `StreamReader`（ストリーム取得・時刻整合） / `parse_data_stream`（バッファ管理・描画）
- 既定サンプリング周波数: `SAMPLING_RATE`（ソース参照: 公開時は 100 Hz）
- バッファ長: `MAX_SAMPLES = duration * SAMPLING_RATE`

## 前提（推奨環境）
- Ubuntu（本ドキュメントは Ubuntu を想定）
- Python 3.12
- 必要ライブラリ: `matplotlib`, `scipy`, `numpy`（バージョンは `pyproject.toml` を参照）

## インストール（Ubuntu）
1. 仮想環境作成・有効化:
```sh
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
```

2. 依存関係のインストール（いずれか）:
- プロジェクトインストール（`pyproject.toml` がある場合、推奨）
```sh
pip install -e .
```
- 最低限の主要ライブラリを直接インストール
```sh
pip install matplotlib scipy numpy
```
- `requirements.txt` を生成したい場合（ローカルで依存ライブラリを入れた後）
```sh
pip freeze > requirements.txt
```

## 実行例（基本形）
shmdump の出力をパイプで渡す基本形:
```sh
shmdump -tq <shm_id> -[L|H|B] <fliter> -f <channel_table> - | python3 shmfall.py -f <channel_table> [オプション]
```

### 必須引数
- `-f, --channel_table` : チャネル定義ファイル（例: `dastest_100ch.tbl`）

### 主なオプション
- `-v, --vminmax` : カラースケール範囲 `vmin:vmax`（例: `-v -20000:20000`）
- `-d, --duration` : 表示期間（秒）。内部で `SAMPLING_RATE` を用いて `MAX_SAMPLES` を決定
- `-n, --normalize` : 正規化モード（`0`=Z-score, `1`=Robust Z-score）
- `-e, --envelope`  : 信号処理（`0`=Envelope(Hilbert), `1`=RMS）
- `-i, --ignore_zeros` : 全ゼロチャネル（欠損）を無視
- `-s, --sort_channels` : チャネルをでソート(`0`=経度, `1`=リファレンスポイントからの距離)
- `-r, --ref_point` : リファレンスポイントの緯度・経度 (緯度/経度)
- `--debug` : 初回フレームを PNG 保存（`shmfall_latest_debug.png`）
- `--help` : 詳細な引数は argparse のヘルプを参照

## 使用例（workspace の example.sh から抜粋）
- 共有メモリ（ID 11）の常時モニタリング例:
```sh
shmdump -tq 11 -H 2.0:1.0:0.5:5.0 -f dastest_100ch.ch - | python3 shmfall.py -i -s0 -f dastest_100ch.ch -v-20000:20000 -d30
```

- エンベロープ処理を有効にしてデバッグ画像出力:
```sh
shmdump -tq 11 -H 2.0:1.0:0.5:5.0 -f dastest_100ch.ch - | python3 shmfall.py -e0 -i -s0 -f dastest_100ch.ch -v0:30000 -d30 --debug
```

- WIN 形式ファイルを shmdump に渡す一例（ファイルをパイプ入力）:
```sh
cat sample_data/25051921.1925 | shmdump -tq -H 2.0:1.0:0.5:5.0 -f dastest_100ch.ch - | python3 shmfall.py -i -s0 -f dastest_100ch.ch -v-20000:20000 -d30 --debug
```

## デバッグとトラブルシュート
- `--debug` を指定すると初回フレームを `shmfall_latest_debug.png` に保存します（描画確認用）。
- チャンネルテーブルの列配置（ID や緯度・経度の列番号）は `shmfall.py` の期待に合わせてください（ソースが split() で解析している列番号に依存します）。
- チャンネルテーブルの作り方はWINシステムのマニュアルページから、[win](https://wwweic.eri.u-tokyo.ac.jp/WIN/man.ja/win.html)の項目の"2.2 パラメタファイルの設定/[2]チャネル表ファイル"を参照してください。
- データが同期しない／ギャップが多い場合は `shmdump` 側のオプションやパイプ経路を確認してください。

## 注意事項・補足
- 本ツールは `shmdump` が出力するテキストフォーマットに依存します。別フォーマットを渡す場合は `shmfall.py` のパーサを調整してください。
- 大規模チャンネル数・高サンプリングでは描画負荷が高くなります。必要に応じて `-d`（duration）やサンプリング設定を調整してください。
- 本文献で用いたチャンネルはこのリポジトリに含まれません。図の出力のみ試したい場合は、任意のチャンネル番号（例えば 0x0001 ~ 0x0064）でチャンネルテーブルを生成してください。
- サンプルデータはいずれも本文中に示す表示例の100chのデータです。

## 参考
- 卜部　卓・束田進也，1992，win─微小地震観測網波形験測支援のためのワークステーション・プログラム（強化版），日本地震学会講演予稿集1992年度秋季大会，P 41.
- 卜部 卓，1994，多チャンネル地震波形データのための共通フォーマットの提案 , 日本地震学会講演予稿集 , No. 2, P24.
- `shmfall.py`（実装本体、StreamReader / parse_data_stream を参照）
- `example.sh`（実行例の収集元）
- `pyproject.toml`（プロジェクト依存の定義）


