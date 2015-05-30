package calculator

func subProfile(a, b []byte) []float64 {
	var subs []float64
	for i := 0; i < len(a); i++ {
		var v float64 = 0
		if a[i] != b[i] {
			v = 1.0
		}
		subs = append(subs, v)
	}

	return subs
}
