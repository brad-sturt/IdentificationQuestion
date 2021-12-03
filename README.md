# IdentificationQuestion

This repository contains the code used in the numerical experiments in the paper:

	B. Sturt (2021). The Value of Robust Assortment Optimization Under Ranking-based Choice Models. Available at SSRN: 


## Citation

If you use the code in this repository in your own research, please cite the above paper as follows:

```
@article{sturt2021identification,
	title={The Value of Robust Assortment Optimization Under Ranking-based Choice Models},
	author={Sturt, Bradley},
	year={2021},
}
```


## Organization

This repository is organized into the following directories: 

* `code/` Contains the scripts for performing the numerical experiments from Sections 4.1 and 4.2. The output of executing the scripts is saved in the data directory. 
	* `revenue_ordered_assortments.jl` Performs all of the numerical experiments with revenue-ordered assortments in Section 4.1. After loading the script into the Julia terminal, the code is executed using the main function. 
	* `two_assortments.jl` Performs the first set of numerical experiments with two assortments from Section 4.2. Specifically, this file generates the data for Figure 5. After loading the script into the Julia terminal, the code is executed using the main function. 
	* `two_assortments_speed_test.jl` Performs the second set of numerical experiments with two assortments from Section 4.2. Specifically, this file generates the data for Figure 6. After loading the script into the Julia terminal, the code is executed using the main function. 
	* `utils.jl` Contains helper functions for the other scripts. 
*  `data/` Contains the output of the experiment scripts. 
	*  `revenue_ordered_assortments_n_equal_4.csv` The output of the revenue_ordered_assortments.jl script.
	*  `two_assortments_n_10_K_10.csv` The output of the two_assortments.jl script. 
	*  `two_assortments_speed.csv` The output of the two_assortments_speed_test.jl script.
* `r_scripts/` Contains the scripts for generating the figures from Sections 4.1 and 4.2.
	* `section_41.R` Contains the code for generating Figures 2 and 3. 
	* `section_42.R` Contains the code for generating Figures 5 and 6.


## License

MIT License

Copyright (c) 2021 Bradley Sturt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


