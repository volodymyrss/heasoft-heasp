An extension of HEASoft heasp module, adding a function to initialize and write an OGIP-compliant spectrum:

```python
import heaspa

heaspa.PHA([0,1],[1,1],10).write("a.fits")
```
