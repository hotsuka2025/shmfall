# shmfall_debugview.py
# デバッグ用：shmfall.py から渡された debug_matrix (time x channel) を
# 任意範囲で切り出し、縦軸を “秒” 単位でラベルして図を保存する。

import os
import numpy as np
import matplotlib.pyplot as plt


# ▼ デフォルトの切り出し範囲（None = 全範囲）
#   y: 時間方向（サンプルインデックス）
#   x: チャンネル方向（インデックス）
#
#   例:
#   DEFAULT_Y_RANGE = (-500, None)  # 最後の500サンプルだけ
#   DEFAULT_X_RANGE = (10, 60)      # ch=10〜60だけ
DEFAULT_Y_RANGE = (500, 1500)
DEFAULT_X_RANGE = (0, 30)


def _make_slice(max_n, range_tuple):
    """range_tuple をスライスに変換する。"""
    if range_tuple is None:
        return slice(None)
    start, stop = range_tuple
    return slice(start, stop)


def plot_matrix_slice(
    debug_matrix: np.ndarray,
    sampling_rate: float,
    y_range=None,
    x_range=None,
    vmin=None,
    vmax=None,
    cbar_label="Amplitude",
    cmap="seismic",
    output_path: str = "shmfall_debug_slice.png",
    show: bool = False,
):
    """
    shmfall.py update_matrix() から渡された最終描画マトリクス debug_matrix を
    任意範囲で切り出し、縦軸を秒単位でラベルした図を出力する。

    Parameters
    ----------
    debug_matrix : np.ndarray
        形状 (n_time, n_chan) の 2D 配列。行方向が時間サンプル。
    sampling_rate : float
        サンプリング周波数 [Hz]。100 Hz なら 1サンプル = 0.01 秒。
    y_range : (start, stop) or None
        行方向（サンプルインデックス）のスライス範囲。
        None の場合は DEFAULT_Y_RANGE を使用。
    x_range : (start, stop) or None
        列方向（チャンネルインデックス）のスライス範囲。
        None の場合は DEFAULT_X_RANGE を使用。
    output_path : str
        出力 PNG ファイル名。
    show : bool
        True なら plt.show() で表示も行う。
    """

    if debug_matrix.ndim != 2:
        raise ValueError("debug_matrix must be 2D array")

    n_time, n_chan = debug_matrix.shape

    if y_range is None:
        y_range = DEFAULT_Y_RANGE
    if x_range is None:
        x_range = DEFAULT_X_RANGE

    y_slice = _make_slice(n_time, y_range)
    x_slice = _make_slice(n_chan, x_range)

    sub = debug_matrix[y_slice, x_slice]
    if sub.size == 0:
        raise ValueError("Selected slice is empty. Check y_range/x_range.")

    # --- サンプルインデックス → 秒への変換 ---
    # スライスの開始・終了インデックス（None の場合は 0 / n_time）
    y_start_idx = 0 if y_slice.start is None else y_slice.start
    y_stop_idx  = n_time if y_slice.stop is None else y_slice.stop

    # 秒へ換算
    y_start_sec = y_start_idx / sampling_rate
    y_stop_sec  = y_stop_idx  / sampling_rate

    # imshow の extent:
    # extent = (xmin, xmax, ymin, ymax)
    # origin='upper' の場合、ymax が上端、ymin が下端として描画される。
    extent = (
        0,                  # x: 左端（チャンネル index 0）
        sub.shape[1],       # x: 右端（チャンネル index max）
        y_stop_sec,         # y: 下端（秒）
        y_start_sec,        # y: 上端（秒）
    )

    # --- プロット ---
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    im = ax.imshow(
        sub,
        aspect="auto",
        interpolation="spline16",
        extent=extent,
        origin="upper",
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )

    ax.set_xlabel("Channel")
    ax.set_ylabel("Elapsed time [sec]")
    fig.colorbar(im, ax=ax, label=cbar_label)

    fig.tight_layout()
    fig.savefig(output_path)
    if show:
        plt.show()
    plt.close(fig)

    print(f"[debugview] saved slice image to {os.path.abspath(output_path)}")
