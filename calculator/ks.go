package calculator

import (
	"github.com/mingzhi/gomath/stat/desc/meanvar"
)

type Ks struct {
	*meanvar.MeanVar
}

func (k *Ks) Append(k2 *Ks) {
	k.MeanVar.Append(k2.MeanVar)
}

func NewKs() *Ks {
	kc := Ks{}
	kc.MeanVar = meanvar.New()
	return &kc
}

func CalcKs(sequences [][]byte) *Ks {
	ks := NewKs()
	for i := 0; i < len(sequences); i++ {
		for j := i + 1; j < len(sequences); j++ {
			subs := Compare(sequences[i], sequences[j])
			for k := 0; k < len(subs); k++ {
				ks.Increment(subs[k])
			}
		}
	}

	return ks
}
