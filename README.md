# Jana2006_Analysis
Python functions for reading refinement data from Jana2006

**Version 0.1**

By Dan Porter, Diamond Light Source
2019

#### How to use:
- Start your refinement in Jana2006
- Find the *.m40 file in the refinement directory
- in a python console:

```python
import Jana2006_Analysis as jana
ref = jana.Refine("refinementfile.m40")
ref.update()
```

- This will print the current refinement results and save them to a text file in your refinement table
- You can also create a latex refinement table:
```python
ref.create_table()
```
