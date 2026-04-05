import logging
from pathlib import Path
from strkit.convert.constants import FORMAT_TRF, FORMAT_BED4
from strkit.convert.converter import convert

TRF_CONVERT_SORT_DAT = Path(__file__).parent / "data" / "trf_convert_sort.dat"

logger = logging.getLogger(__name__)


def test_convert_trf_dat_to_bed4(capsys):
    res = convert(str(TRF_CONVERT_SORT_DAT), FORMAT_TRF, FORMAT_BED4, sort=False, logger=logger)
    assert res == 0
    captured = capsys.readouterr()
    assert captured.out == "Chr_Y.H1\t569\t609\tGGGCGGCGGGGGAAAAAAG\nChr_Y.H1\t0\t62\tTATCACCGGGGAGGAGGAGGGTAGTCAG\n"


def test_convert_trf_dat_to_bed4_sorted(capsys):
    res = convert(str(TRF_CONVERT_SORT_DAT), FORMAT_TRF, FORMAT_BED4, sort=True, logger=logger)
    assert res == 0
    captured = capsys.readouterr()
    assert captured.out == "Chr_Y.H1\t0\t62\tTATCACCGGGGAGGAGGAGGGTAGTCAG\nChr_Y.H1\t569\t609\tGGGCGGCGGGGGAAAAAAG\n"
