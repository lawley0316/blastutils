# blastuitls：BLAST工具集

在使用BLAST时，经常需要对比对结果进行操作，例如：读取、排序、过滤等。 Biopython是一个不错的选择，
但是新版Biopython貌似要放弃对tab分隔的文本文件的支持，进而只支持XML、XML2格式，这会带来一些问题：

（1）XML文本可读性很差；
（2）XML文件比tab分隔的文本文件大很多，我碰到过一次比对结果有1.5TB，而tab分隔的文本文件只有40GB；
（3）`Bio.Blast.parse`有bug，不知道啥时候就抛出`CorruptedXMLError`（**这点很关键**）；
（4）biopython对排序、过滤等操作的支持有点弱。

## 安装

```shell
pip install blastutils
```

## `HSP`、`Hit`、`Record`

### `HSP`

在BLAST中，HSP指的是两个序列之间的一段局部比对区域，这段比对区域在某些特定的评估标准（如得分、e-value、序列覆盖度等）上具有较高的相似性。
HSP通常是BLAST比对中最重要的部分，因为它们代表了两个序列之间最强的相似区域。HSP对象有11个属性:

1. `qstart`，比对在查询序列上的起始位置；
2. `qend`，比对在查询训练上的终止位置；
3. `sstart`，比对在目标序列上的起始位置；
4. `send`，比对在目标训练上的终止位置；
5. `mismatch`，错配数；
6. `gapopen`，gap数；
7. `length`，比对长度；
8. `pident`，比对的identity；
9. `qcovhsp`，<float>, HSP在查询序列上的覆盖度；
10. `bitscore`，<float>, Bit score；
11. `evalue`，<float>, Expect value。

### `Hit`

在BLAST中，Hit是指查询序列与数据库中某个目标序列之间的比对结果，通常包含一个或多个HSP。Hit对象有3个属性：

1. `sseqid`，目标序列编号；
2. `slen`，目标序列长度；
3. `hsps`，相应的HSPs。

### `Record`

一个Record代表一个查询序列的所有比对结果，通常包含一个或多个Hit。Record对象有3个属性：

1. `qseqid`，查询序列编号；
2. `qlen`，查询序列长度；
3. `hits`，相应的Hits。

## 使用

### 读写BLAST结果

注意：使用`blastutils`读取BLAST结果需要将比对时`outfmt`参数的值设置为：
`6 qseqid qlen sseqid slen qstart qend sstart send mismatch gapopen length pident qcovhsp bitscore evalue`

```python
from blastutils import BlastOutputFile, Reader, Writer

with BlastOutputFile('example.txt') as file:
    reader = Reader(file)
    records = []
    for record in reader:
        records.append(record)

with BlastOutputFile('new-example.txt', 'w') as file:
    writer = Writer(file)
    for record in records:
        writer.write(record)
```

`BlastOutputFile`是对BLAST结果文件的封装，有两个参数：

（1）`path`，文件路径；
（2）`mode`，打开模式，`r`代表读取（默认），`w`代表写入。

可以使用像上面使用上下文管理使用也可以手动管理：

```python
from blastutils import BlastOutputFile, Reader

file = BlastOutputFile('example.txt')
file.open()
row = file.read()  # ['seq1', 100, 'ref1', 100, 1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 8.22e-50]
file.close()
```

把`BlastOutputFile`对象作为参数传入`Reader`中即可进行读取。`Reader`对象是一个迭代器，可以像上面使用`for`循环进行迭代，
也可以使用`next`函数进行操作：

```python
from blastutils import BlastOutputFile, Reader

file = BlastOutputFile('example.txt')
file.open()
reader = Reader(file)
record1 = next(reader)
record2 = next(reader)
file.close()
```

### 过滤

过滤是一个非常常见的操作，例如过滤掉相似度过低、覆盖度过低的Hit是很多生物信息流程一上来就要做的。但是过滤的规则却是多种多样的，
不可能设置一套普适的过滤方法。所以`blastutils`并没有内置过多的过滤规则，更多的是提供API接口让用户自己去实现适合的过滤规则。

`blastutils`定义了一个`Filter`基类，其所有子类都都必须实现`__call__`方法，该方法传入一个Hit对象，如果返回`True`代表该Hit需要被保留，否则就删除。

例如，可以这样实现一个过滤规则，当该Hit的最好的那个HSP（对于HSP，BLAST的输出结果是排好序的，第一个就是最好的）的pident小于90的都不要，就可以这样实现：

```python
from blastutils import Filter

class MinSimilarity(Filter):
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, hit):
        if hit.is_empty():  # 对于没有HSP的Hit肯定要过滤掉
            return False
        return hit.hsps[0].pident >= self.threshold

record.filter(MinHSPSimilarity(90))  # 过滤掉小于90的Hit
```

### 排序

当获取到BLAST结果后，一个最常见的需求就是想知道每个查询序列的最佳匹配是哪个，常见的做法就是对所有Hits进行排序，然后取第一个。
同过滤一样，排序的规则同样是多种多样的，也没有一套排序规则是普适的。

`blastutils`定义了一个`Compare`基类，其所有子类都必须实现`__call__`方法，该方法传入两个Hit对象hit1和hit2：

（1）返回-1（或任意一个小于0的数）代表hit1好于hit2；
（2）返回0代表hit1与hit2同样好；
（3）返回1（或任意一个大于0的数）代表hit2好于hit1。

例如，可以实现这样一个排序规则：

（1）evalue约小越好；
（2）evalue相同时，bitscore越大越好。

```python
from blastutils import Compare

class ByEvalueBitscore(Compare):
    def __call__(self, hit1, hit2):
        if hit2.is_empty():  # hit2没有HSP，肯定hit1好
            return -1
        if hit1.is_empty():  # hit1没有HSP，肯定hit2好
            return 1
        hsp1 = hit1.hsps[0]
        hsp2 = hit2.hsps[0]
        if hsp1.evalue < hsp2.evalue:
            return -1
        if hsp1.evalue > hsp2.evalue:
            return 1
        return -1 if hsp1.bitscore >= hsp2.bitscore else 1

record.best(ByEvalueBitScore())  # 没有排好序，依然可以获取最佳Hit
record.sort(ByEvalueBitScore())  # 排序
```