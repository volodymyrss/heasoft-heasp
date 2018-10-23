An extension of HEASoft heasp module, adding a function to initialize and write an OGIP-compliant spectrum:

```python
import heaspa

heaspa.PHA(counts=[0,1],count_errors=[1,1],exposure=10).write("a.fits")
```
