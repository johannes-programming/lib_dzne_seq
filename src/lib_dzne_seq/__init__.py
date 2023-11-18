import Bio.Seq as _Seq
import Bio.SeqIO as _SeqIO
import Bio.SeqRecord as _SeqRecord
import lib_dzne_filedata as _fd
import lib_dzne_math.na as _na


class SeqRead(_fd.FileData):
    def __len__(self):
        return len(self.seq)
    def __iter__(self):
        return iter(self.seq)
    def __getitem__(self, key):
        digest = self._digest_data(self.data)
        digest['seq'] = digest['seq'][key]
        return self._calc_data(**digest)
    def __setitem__(self, key, value):
        digest = self._digest_data(value)
        digest['seq'][key] = value
        self.data = self._calc_data(**digest)
    def __delitem__(self, key):
        digest = self._digest_data(self.data)
        del digest['seq'][key]
        self.data = self._calc_data(**digest)
    @classmethod
    def _load(cls, /, file):
        return _SeqIO.read(file, cls._format)
    def _save(self, /, file):
        _SeqIO.write(file, type(self)._format, self.data)
    @classmethod
    def _default(cls):
        return cls._calc_data(seq="", qv=0)
    @classmethod
    def from_file(cls, file, /):
        parentclass, = cls.__bases__
        return parentclass.from_file(file, PHDRead, ABIRead)
    def to_TOMLData(self):
        return _fd.TOMLData(self._digest_data(self.data))
    @classmethod
    def clone_data(cls, data, /):
        return cls._calc_data(**cls._digest_data(data))
    @classmethod
    def _calc_data(cls, *, qv, seq):
        qv = int(round(qv))
        if (qv < 0) or (qv > 100):
            raise ValueError()
        seq = normstr(seq)
        seq = _Seq.Seq(seq)
        ans = _SeqRecord.SeqRecord(seq)
        ans.letter_annotations['phred_quality'] = [qv] * len(ans)
        return ans
    @classmethod
    def _digest_data(cls, data):
        seq = normstr(data.seq)
        seq = _Seq.Seq(seq)
        ll = data.letter_annotations['phred_quality']
        if len(ll) == 0:
            return dict(seq=seq, qv=0)
        qv = int(round(sum(ll) / len(ll)))
        if (qv < 0) or (qv > 100):
            raise ValueError()
        return dict(seq=seq, qv=qv)
    @property
    def seq(self):
        kwargs = self._digest_data(self.data)
        return kwargs['seq']
    @seq.setter
    def seq(self, value):
        kwargs = self._digest_data(self.data)
        kwargs['seq'] = value
        self.data = self._calc_data(**kwargs)
    @property
    def qv(self):
        kwargs = self._digest_data(self.data)
        return kwargs['qv']
    @qv.setter
    def qv(self, value):
        kwargs = self._digest_data(self.data)
        kwargs['qv'] = value
        self.data = self._calc_data(**kwargs)
    @property
    def tr(self):
        return tr(self.seq)
    @property
    def contains_stop(self):
        return '*' in self.tr
    @property
    def gravy(self):
        return gravy(self.tr)

class PHDRead(SeqRead):
    _ext = '.phd'
    _format = 'phd'
class ABIRead(SeqRead):
    _ext = '.ab1'
    _format = 'abi'


class FASTAData(_fd.FileData):
    _ext = '.fasta'
    def __len__(self):
        return len(self.data)
    def __iter__(self):
        return iter(self.data)
    def __add__(self, other):
        return self._add(self, other)
    def __radd__(self, other):
        return self._add(other, self)
    def __mul__(self, other):
        return self._mul(self, other)
    def __rmul__(self, other):
        return self._mul(self, other)
    def __getitem__(self, key):
        ans = self.data[key]
        if type(ans) is list:
            ans = type(self)(ans)
        return ans
    def __setitem__(self, key, value):
        data = self.data
        data[key] = value
        self.data = data
    def __delitem__(self, key):
        data = self.data
        del data[key]
        self.data = data
    @classmethod
    def _add(cls, *objs):
        items = list()
        for obj in objs:
            dataObj = cls(obj)
            items += dataObj.data
        ans = cls(items)
        return ans
    @classmethod
    def _mul(cls, dataObj, n):
        return cls(dataObj.data * n)
    @classmethod
    def clone_data(cls, data):
        ans = list()
        for item in data:
            ans.append(cls.clone_item(item))
        return ans
    @classmethod
    def _load(cls, file):
        return _SeqIO.parse(
            handle=file, 
            format='fasta',
        )
    def _save(self, file):
        _SeqIO.write(
            handle=file, 
            format='fasta', 
            sequences=self.data,
        )
    @staticmethod
    def _digest_item(item):
        ans = dict()
        ans['seq'] = _Seq.Seq(normstr(item.seq))
        ans['name'] = str(item.id)
        return ans
    @staticmethod
    def _calc_item(*, name, seq):
        ans = _SeqRecord.SeqRecord(
            seq=_Seq.Seq(normstr(seq)),
            id=str(name),
            description="",
        )
        return ans
    @classmethod
    def clone_item(cls, item):
        digest = cls._digest_item(item)
        ans = cls._calc_item(**digest)
        return ans
    @staticmethod
    def _default():
        return list()
    def ids(self):
        return [rec.id for rec in self.data]
    def seqs(self):
        return [rec.seq for rec in self.data]
    def append(self, name, seq):
        data = self.data
        item = (name, seq)
        data.append(item)
        self.data = data
    def pop(self, *args, **kwargs):
        data = self.data
        ans = data.pop(*args, **kwargs)
        self.data = data
        return ans
    def clear(self):
        self.data = list()




def data(*, seq, go, end):
    if _na.anyisna(seq, go, end):
        return None
    ans = dict()
    ans['go'] = go
    ans['end'] = end
    ans['seq'] = cut(seq, go, end, strict=True)
    ans['seq_len'] = len(ans['seq'])
    ans['tr'] = tr(ans['seq'])
    ans['tr_len'] = len(ans['tr'])
    ans['contains_stop'] = '*' in ans['tr']
    ans['gravy'] = gravy(ans['tr'])
    return ans

def seq3(seq, go=None, end=None):
    seq = _Seq.Seq(seq)
    seq = cut(seq, go=go, end=end)
    while len(seq) % 3:
        seq += 'N'
    return seq

def cut(seq, go=None, end=None, strict=False):
    if go is not None and end is not None and go > end:
        raise IndexError(f"One cannot cut the seq {seq.__repr__()} from {go} until {end}! ")
    if go is not None and go < 0:
        if strict:
            raise IndexError(f"The value go={go} will not be accepted! ")
        prefix = 'N' * (0 - go)
        go = None
    else:
        prefix = ""
    if end is not None and end > len(seq):
        if strict:
            raise IndexError(f"The value end={end} will not be accepted (len={len(seq)})! ")
        suffix = 'N' * (len(seq) - end)
        end = None
    else:
        suffix = ""
    if None not in {go, end}:
        if go > end:
            raise IndexError(f"The values go={go} and end={end} are incompatible! ")
    return prefix + seq[go:end] + suffix

def tr(seq, go=None, end=None):
    return str(seq3(seq, go=go, end=end).translate())


def normstr(seq):
    seq = str(seq).upper()
    forbidden = set(seq) - set("ACGTUN-")
    if len(forbidden):
        raise ValueError(
            f"The letters {forbidden} are forbidden! "
        )
    return seq


def gravy(translation):
    values = {
        'A':1.8,
        'C':2.5,
        'D':-3.5,
        'E':-3.5,
        'F':2.8,
        'G':-0.4,
        'H':-3.2,
        'I':4.5,
        'K':-3.9,
        'L':3.8,
        'M':1.9,
        'N':-3.5,
        'P':-1.6,
        'Q':-3.5,
        'R':-4.5,
        'S':-0.8,
        'T':-0.7,
        'V':4.2,
        'W':-0.9,
        'X':float('nan'),
        'Y':-1.3,
        '-':float('nan'),
        '*':float('nan'),
    }
    answers = [values[k] for k in translation if _na.notna(values[k])]
    if len(answers):
        return sum(answers) / len(answers)
    return float('nan')


    




 
