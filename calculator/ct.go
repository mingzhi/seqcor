package calculator

import (
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc"
)

type AutoCovFFTW struct {
	N            int
	maskXYs, xys []float64
	mean         *desc.Mean
	dft          *correlation.FFTW
}

func NewAutoCovFFTW(n int, dft *correlation.FFTW) *AutoCovFFTW {
	var c AutoCovFFTW
	c.dft = dft
	c.N = n
	c.maskXYs = make([]float64, c.N)
	c.xys = make([]float64, c.N)
	c.mean = desc.NewMean()

	return &c
}

func (c *AutoCovFFTW) Increment(xs []float64) {
	masks := make([]float64, len(xs))
	for i := 0; i < len(masks); i++ {
		masks[i] = 1.0
		c.mean.Increment(xs[i])
	}

	maskXYs := c.dft.AutoCorr(masks)
	xys := c.dft.AutoCorr(xs)

	for i := 0; i < len(c.xys); i++ {
		c.xys[i] += (xys[i] + xys[(len(xys)-i)%len(xys)])
		c.maskXYs[i] += (maskXYs[i] + maskXYs[(len(maskXYs)-i)%len(maskXYs)])
	}
}

func (c *AutoCovFFTW) Append(c2 *AutoCovFFTW) {
	c.mean.Append(c2.mean)
	for i := 0; i < len(c2.xys); i++ {
		c.xys[i] += c2.xys[i]
		c.maskXYs[i] += c2.maskXYs[i]
	}
}

func (c *AutoCovFFTW) GetResult(i int) float64 {
	pxy := c.xys[i] / c.maskXYs[i]
	pxpy := c.mean.GetResult() * c.mean.GetResult()
	return pxy - pxpy
}

type AutoCovFFT struct {
	N            int
	maskXYs, xys []float64
	mean         *desc.Mean
	circular     bool
}

func NewAutoCovFFT(maxl int, circular bool) *AutoCovFFT {
	var c AutoCovFFT
	c.N = maxl
	c.maskXYs = make([]float64, c.N)
	c.xys = make([]float64, c.N)
	c.mean = desc.NewMean()
	c.circular = circular

	return &c
}

func (c *AutoCovFFT) Increment(xs []float64) {
	masks := make([]float64, len(xs))
	for i := 0; i < len(masks); i++ {
		masks[i] = 1.0
		c.mean.Increment(xs[i])
	}

	maskXYs := correlation.AutoCorrFFT(masks, c.circular)
	xys := correlation.AutoCorrFFT(xs, c.circular)

	for i := 0; i < len(c.xys); i++ {
		c.xys[i] += (xys[i] + xys[(len(xys)-i)%len(xys)])
		c.maskXYs[i] += (maskXYs[i] + maskXYs[(len(maskXYs)-i)%len(maskXYs)])
	}
}

func (c *AutoCovFFT) Append(c2 *AutoCovFFT) {
	c.mean.Append(c2.mean)
	for i := 0; i < len(c2.xys); i++ {
		c.xys[i] += c2.xys[i]
		c.maskXYs[i] += c2.maskXYs[i]
	}
}

func (c *AutoCovFFT) GetResult(i int) float64 {
	pxy := c.xys[i] / c.maskXYs[i]
	pxpy := c.mean.GetResult() * c.mean.GetResult()
	return pxy - pxpy
}

type AutoCov struct {
	N     int
	corrs []*correlation.BivariateCovariance
}

func NewAutoCov(maxl int, bias bool) *AutoCov {
	cc := AutoCov{}
	cc.corrs = make([]*correlation.BivariateCovariance, maxl)
	for i := 0; i < maxl; i++ {
		cc.corrs[i] = correlation.NewBivariateCovariance(bias)
	}
	cc.N = maxl
	return &cc
}

func (cc *AutoCov) Increment(i int, x, y float64) {
	cc.corrs[i].Increment(x, y)
}

func (cc *AutoCov) GetResult(i int) float64 {
	return cc.corrs[i].GetResult()
}

func (cc *AutoCov) GetMeanXY(i int) float64 {
	return cc.corrs[i].MeanX() * cc.corrs[i].MeanY()
}

func (cc *AutoCov) GetN(i int) int {
	return cc.corrs[i].GetN()
}

func (cc *AutoCov) Append(cc2 *AutoCov) {
	for i := 0; i < len(cc.corrs); i++ {
		cc.corrs[i].Append(cc2.corrs[i])
	}
}

func CalcCtFFTW(sequences [][]byte, dft *correlation.FFTW) *AutoCovFFTW {
	n := len(sequences[0])
	ct := NewAutoCovFFTW(n, dft)
	for i := 0; i < len(sequences); i++ {
		for j := i + 1; j < len(sequences); j++ {
			subs := Compare(sequences[i], sequences[j])
			ct.Increment(subs)
		}
	}

	return ct
}

func CalcCtFFT(sequences [][]byte, maxl int, circular bool) *AutoCovFFT {
	ct := NewAutoCovFFT(maxl, circular)
	for i := 0; i < len(sequences); i++ {
		for j := i + 1; j < len(sequences); j++ {
			subs := Compare(sequences[i], sequences[j])
			ct.Increment(subs)
		}
	}

	return ct
}

func CalcCt(sequences [][]byte, maxl int) *AutoCov {
	ct := NewAutoCov(maxl, false)
	for i := 0; i < len(sequences); i++ {
		for j := i + 1; j < len(sequences); j++ {
			subs := Compare(sequences[i], sequences[j])

			for l := 0; l < maxl; l++ {
				for k := 0; k < len(subs)-l; k++ {
					x, y := subs[k], subs[k+l]
					ct.Increment(l, x, y)
				}
			}
		}
	}

	return ct
}
