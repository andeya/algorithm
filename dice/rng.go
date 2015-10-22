package dice

import (
	"sync"
	"math/rand"
	"time"
)

type lockedSource struct {
	lk  sync.Mutex
	src rand.Source
}

func (r *lockedSource) Int63() (n int64) {
	r.lk.Lock()
	n = r.src.Int63()
	r.lk.Unlock()
	return
}

func (r *lockedSource) Seed(seed int64) {
	r.lk.Lock()
	r.src.Seed(seed)
	r.lk.Unlock()
}

var localRng = rand.New(&lockedSource{src: rand.NewSource(time.Now().UnixNano())})

type intRng interface {
	Intn(n int) int
}
