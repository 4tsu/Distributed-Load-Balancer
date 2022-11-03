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
[6]:https://www.lammps.org/



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
md->set_config("ファイル名");
```
book keeping法については[こちら][5]を参照。

粒子配置の読込ファイルについては、LAMMPSの.dumpファイルのフォーマットに準拠し、`x, y, z, vx, vy, vz`を入力する形となる。`sample/sample.dump`がお手本となる。また、3次元のシミュレーション結果から、断面を切り出して2次元粒子配置として読み込ませるコードが`sample/ttt.py`である。このあたり、詳しくは[サンプル](#サンプルsample)を参照。

一方で、ファイル名を指定しない場合は、`set_config("make");`とすることで、粒子を一様に配置した系を生成する。

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



## サンプル（sample）
粒子配置ファイルを生成するサンプルが、sampleフォルダ内にある。分子動力学シミュレーションソフト[LAMMPS][6]の出力、`.dump`形式を粒子配置入力の手本にしているため、LAMMPSを途中まで走らせて、続きから本プログラムで走らせる事が可能（本プログラムの性能評価のときに便利だった）。

- sample.dump : `.dump`ファイルのサンプル。本プログラムは、この形式の粒子配置ファイルを読み込むことができる。自分で記述するときは、`ITEM: TIMESTEP`の値を0にしてはいけない。ここが0だと、自動で飛ばすように組んでいる（LAMMPSの続きからやるときに便利）。`ITEM: TIMESTEP`は書かなくても読み込まれ、その場合は、step 0 からシミュレーションがスタートする。
- run.in : LAMMPSを使うスクリプト。落ち着いた状態の粒子配置と速度の情報を得るために、初期緩和だけLAMMPSを使用できる。本プログラムは温度制御を未実装なので、所望の温度に落ち着くまでLAMMPSを使用したりする。
- makeconf.py : 液滴や、気泡の原子配置を作り、LAMMPSに与えることができる。これは本プログラムでは読み込めない。
- ttt.py : 3次元シミュレーションの結果得た`.dump`ファイルから一部の粒子を抜き出し、2次元の粒子配置ファイルを生成する。

`makeconf.py` -> `run.in` -> `ttt.py` の順に実行すると、デフォルトでは`smpl2d.dump`がsampleフォルダ内に生成されるので、これを`md.exe`と同じ場所に移動し、`set_config()`の引数にファイル名を指定することで、液滴のサンプルシミュレーションが動く。
```
cd sample
python3 makeconf.py

mpirun -np 並列数 lmp_mpi -in run.in
# or 
lmp_serial -in run.in

python3 ttt.py

mv smpl2d.dump ..

vi src/main.cpp
# md->set_config("smpl2d.dump");に書き換える

mpirun -np 並列数 ./md.exe
```

## バージョン
mainブランチ参照。

または、[各バージョンの性能](performance.md)を参照。
