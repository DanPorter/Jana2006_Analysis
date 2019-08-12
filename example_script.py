"""
Example script for Jana2006_Analysis
"""

import Jana2006_Analysis as jana

f = r"C:\Users\dgpor\Dropbox\Mini Projects\Ca2RuO4\XRD Ca2RuO4 Feb 2017\Refine 22 March 17\Pbca\90\Ca2RuO4_Feb_Temp_01_90_Dan.cif"

ref = jana.Refine(f)

#ref.update('A test refinement')
print(ref.refine_results())
#ref.create_table()
