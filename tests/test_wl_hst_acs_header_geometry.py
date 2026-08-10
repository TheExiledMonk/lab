from pathlib import Path
from astropy.io import fits
from pbuf.wl.hst_acs_header_geometry import read_headers

def test_header_reader_never_requests_data(tmp_path,monkeypatch):
    p=tmp_path/"x.fits"
    a,b=fits.ImageHDU(name="SCI"),fits.ImageHDU(name="SCI");a.header["EXTVER"]=1;b.header["EXTVER"]=2
    fits.HDUList([fits.PrimaryHDU(),a,b]).writeto(p)
    monkeypatch.setattr(fits,"getdata",lambda *a,**k: (_ for _ in ()).throw(AssertionError("pixel read")))
    primary,sci=read_headers(p);assert len(sci)==2 and primary is not None
