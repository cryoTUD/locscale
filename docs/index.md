# Overview

<style>
.c-compare {
  position: relative;
  width: 100%;
  max-width: 600px;
  aspect-ratio: 1 / 1;
  overflow: hidden;
}

.c-compare__left,
.c-compare__right {
  position: absolute;
  top: 0;
  left: 0;
  height: 100%;
  width: 100%;
  object-fit: contain;
}

.c-compare__right {
  clip-path: inset(0 0 0 var(--value, 50%));
}

.c-compare__range {
  position: absolute;
  bottom: 10px;
  left: 10px;
  width: 90%;
  z-index: 10;
}
</style>

<div class="c-compare" style="--value:50%;">
  <img class="c-compare__left" src="imgs/emd19995.png"
       alt="Raw map" />
  <img class="c-compare__right" src="imgs/emd19995_fem.png"
       alt="Feature-enhanced map" />
  <input type="range" class="c-compare__range" min="0" max="100" value="50"
         oninput="this.parentNode.style.setProperty('--value', this.value + '%')" />
</div>


## What is LocScale2

The primary purpose ofdlkkfklskflskfkle
<br><br>
What does the 'local scale' 

