import Bio.Seq as _Seq
import Bio.SeqIO as _SeqIO
import Bio.SeqRecord as _SeqRecord
import lib_dzne_filedata as _fd
import lib_dzne_math.na as _na


class SeqRead(_fd.FileData):
    def __len__(self):
        return len(self.seq)
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
        return parentclass.from_file(file, PHDSeqReadData, ABISeqReadData)
    def to_TOMLData(self):
        return _fd.TOMLData(cls._digest_data(self.data))
    @classmethod
    def clone_data(cls, data, /):
        return cls._calc_data(**cls._digest_data(data))
    @classmethod
    def _calc_data(cls, *, qv, seq):
        qv = int(round(qv))
        if (qv < 0) or (qv > 100):
            raise ValueError()
        seq = normstr(seq)
        ans = _SeqRecord.SeqRecord(seq)
        ans.letter_annotations['phred_quality'] = [qv] * len(ans)
        return ans
    @classmethod
    def _digest_data(cls, data):
        seq = normstr(data.seq)
        ll = data.letter_annotations['phred_quality']
        if len(ll) == 0:
            return dict(seq=seq, qv=0)
        qv = int(round(sum(ll) / len(ll)))
        if (qv < 0) or (qv > 100):
            raise ValueError()
        return dict(seq=seq, qv=qv)
    @property
    def seq(self):
        return self.data.seq
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

class PHDRead(SeqReadData):
    _ext = '.phd'
    _format = 'phd'
class ABIRead(SeqReadData):
    _ext = '.ab1'
    _format = 'abi'

def data(*, seq, go, end):
    if _na.anyisna(seq, go, end):
        return None
    ans = dict()
    ans['go'] = go
    ans['end'] = end
    ans['seq'] = cut(seq, go, end, strict=True)
    ans['seq-len'] = len(ans['seq'])
    ans['tr'] = tr(ans['seq'])
    ans['tr-len'] = len(ans['tr'])
    ans['contains-stop'] = '*' in ans['tr']
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
        raise IndexError(f"One cannot cut the seq {ascii(str(seq))} from {go} until {end}! ")
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
    }
    answers = [values[k] for k in translation if _na.notna(values[k])]
    if len(answers):
        return sum(answers) / len(answers)
    return float('nan')


    




 
