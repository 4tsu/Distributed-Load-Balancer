# myMD
C++とMPIで書いた分子動力学シミュレーション（Molecular Dynamics, MD）コードです。自分に必要な機能のみ実装しています。

プロセス並列（分散メモリ）。ノートPCでもスーパーコンピュータでも。



## 参考（Reference）

- [分子動力学法ステップ・バイ・ステップ][1]
- [一週間でなれる！スパコンプログラマ][2]
- 他

[1]:https://github.com/kaityo256/mdstep
[2]:https://github.com/kaityo256/sevendayshpc
[3]:https://github.com/kaityo256/sevendayshpc/blob/main/day1/README.md
[4]:https://polymer.apphy.u-fukui.ac.jp/~koishi/cdview.php
[5]:https://qiita.com/kaityo256/items/cfacbf6f1136de63bd97



## 環境（Environment）

- MPI（[一週間スパコンのDay1][3]を参照）
- C++（C++17以上が好ましい）

C++の`filesystem`が使えるといろいろと便利。なくても動く。



## 使い方（Usage）
### ビルド
一番上の階層で`make`、その場に「`md.exe`」が生成される。

「`filesystem`が使えない」という旨のエラーが出た場合は、[後述のパターン](#filesystemが使えないとき)。

また、環境によっては`mpic++`ではなく`mpicxx`が良い。`makefile`の`CC`を書き換えておく。
```
CC = mpicxx
```

### パラメータの変更
`src/main.cpp`にてシミュレーションパラメータの指定。変更の都度、ビルド。
```
md->set_params([ステップ数], [粒子位置ダンプ間隔(step)], [時間刻み]);
md->set_box([粒子数], [シミュレーションボックスのx長さ], [シミュレーションボックスのy長さ], [カットオフ距離]);
md->set_margin([book keeping法のマージン]);
```
book keeping法については[こちら][5]を参照。

### 実行
```
mpirun -np 並列数 ./md.exe
```
ハイパースレッディングに対応するCPUは、`--oversubscribe`オプションの指定により、並列数の限界が押し上がる。
```
mpirun --oversubscribe -np 並列数 ./md.exe
```
この辺りはMPIの話なので、詳しくは[他の記事][2]などを参照。

### 結果の取得
主に得られる結果
- 各ステップにおける粒子位置のファイル（cdv/.cdv）
- 運動エネルギーとポテンシャルエネルギー、及び系の全エネルギーの変化（energy.dat）

### 結果の可視化
- .cdvは[cdview][4]という可視化プログラムのフォーマットで出力。cdviewを使えば粒子の動きをインタラクティブに見ることができる。また、cdviewが使用できない環境向けに、`vis`フォルダ内の`vis.py`の実行で、`vis`にgifアニメーションが生成される（cdview.gif）。
- energy.datはgnuplot向け。`vis`フォルダ内の`energy.plt`の実行でグラフを得る（energy.png）。
- gnuplotスクリプトは文字コード絡みのエラーを起こしやすい。そこで、`vis.py`で代替できるようにしてある。`vis.py`内で`plot_energy()`関数を実行することで、energy.pngを取得。
```
# vis.py
plot_energy("../energy.dat")
```
- `python3 vis/vis.py`は`make fig`によって実行可能。



### `filesystem`が使えないとき
#### ビルド
`makefile`中の11行目くらい
```
OPTIONS += -DFS
```
をコメントアウトし、`make`。
#### 実行
出力ファイルが前回のシミュレーション結果への追記になってしまう。実行前にこれらを削除しておく。
```
*.cdv
energy.dat
mpirun -np 並列数 ./md.exe
```
#### 結果
得られる結果は[上記](#結果の取得)と同じ。しかし、`cdv`フォルダが作成されず、.cdvファイルは`md.exe`が存在する場所に生成される。`vis.py`による可視化は同様に行うことができる。

## バージョン
mainブランチ参照。
