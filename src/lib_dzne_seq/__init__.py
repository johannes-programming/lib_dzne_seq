import Bio.Seq as _Seq
import Bio.SeqIO as _SeqIO
import Bio.SeqRecord as _SeqRecord
import gravy as _gravy
import lib_dzne_filedata as _fd
import na_quantors as _na
import seqpad as _sp


class SeqRead:
    def __init__(self, **kwargs):
        keys = set(kwargs.keys())
        if keys == {'seq', 'qv'}:
            self.seq = kwargs['seq']
            self.qv = kwargs['qv']
        elif keys == {'record'}:
            self.record = kwargs['record']
        elif keys == {'read'}:
            self.seq = kwargs['read'].seq
            self.qv = kwargs['read'].qv
        elif keys == {'file', 'format'}:
            self.record = _SeqIO.read(kwargs['file'], kwargs['format'])
        elif keys == {}:
            self.seq = ""
            self.qv = 0
        else:
            raise NotImplementedError()
    def __len__(self):
        return len(self.seq)
    def save(self, file, format):
        _SeqIO.write(file, format, self.record)
    @property
    def record(self):
        ans = _SeqRecord.SeqRecord(self.seq)
        ans.letter_annotations['phred_quality'] = [self.qv] * len(self)
        return ans
    @record.setter
    def record(self, value):
        self.seq = value.seq
        ll = value.letter_annotations['phred_quality']
        if len(ll):
            self.qv = 100 * sum(ll) / len(ll)
        else:
            self.qv = 0
    @property
    def seq(self):
        return self._seq
    @seq.setter
    def seq(self, value):
        self._seq = _Seq.Seq(value)
    @property
    def qv(self):
        return self._qv
    @qv.setter
    def qv(self, value):
        value = int(round(value))
        if 0 <= value <= 100:
            self._qv = value
        else:
            raise ValueError("Quality Value must be in the range from 0 to 100. ")
    

class SeqReadData(_fd.FileData):
    @staticmethod
    def clone_data(data):
        return SeqRead(read=data)
    @classmethod
    def _load(cls, /, file):
        return SeqRead(file=file, format=cls._format)
    def _save(self, /, file):
        self._data.save(file=file, format=type(self)._format)
    @classmethod
    def _default(cls):
        return SeqRead()
    @classmethod
    def from_file(cls, file, /):
        return super().from_file(file, PHDSeqReadData, ABISeqReadData)

class PHDSeqReadData(SeqReadData):
    _ext = '.phd'
    _format = 'phd'
class ABISeqReadData(SeqReadData):
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
    seq = _sp.seqpad(seq)
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
    return _gravy.calculate(translation)


    




 
