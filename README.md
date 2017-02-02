# Eclipsing-binaries-library
### Eclisping binaries light curves generator

#### Dependencies:
* python        3
* CGAL          4.2
* numpy      >= 1.11.2
* scipy      >= 0.18.1
* matplotlib >= 1.5.3
* pymysql    >= 0.7.9

#### Getting Started
1. Install dependencies
2. Create MySQL database __elisa_assets__, download and load [database](https://mega.nz/#!7YNUnLBI!d4zzq3HJBPn3LK7H8kHSryPN-xexa6Oh5yVNxbi3WhU)
3. Change default variables (host, user, password) in *__/globe/variables.py__*, default:
   * HOST = "localhost"
   * USER = "root"
   * PWD = "toor"
4. Compile (CGAL 4.2 is necessary) *__/cpp/cgal_build/poisson_reconstruction.cpp__* (use a cmake) and copy __poisson_reconstruction__ to *__/bin/__*


In case of any questions, coments or concerns, please contact me on stealth.on.github[gmail.com].
            
