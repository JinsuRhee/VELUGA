gfortran -fopenmp -fPIC -c find_domain.f90 -o find_domain.o
gcc -fopenmp -fPIC -c find_domain_c2f.c -o find_domain_c2f.o
gcc -shared -fopenmp find_domain.o find_domain_c2f.o -o find_domain.so -lgfortran

gfortran -fopenmp -fPIC -c jsrd_part_totnum.f90 -o jsrd_part_totnum.o
gcc -fopenmp -fPIC -c jsrd_part_totnum_c2f.c -o jsrd_part_totnum_c2f.o
gcc -shared -fopenmp jsrd_part_totnum.o jsrd_part_totnum_c2f.o -o jsrd_part_totnum.so -lgfortran

gfortran -fopenmp -fPIC -c jsrd_part.f90 -o jsrd_part.o
gcc -fopenmp -fPIC -c jsrd_part_c2f.c -o jsrd_part_c2f.o
gcc -shared -fopenmp jsrd_part.o jsrd_part_c2f.o -o jsrd_part.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c js_kdtree.f90 -o js_kdtree.o
gfortran -fopenmp -fPIC -mcmodel=large -c get_contam.f90 -o get_contam.o
gcc -fopenmp -fPIC -c get_contam_c2f.c -o get_contam_c2f.o
gcc -shared -fopenmp get_contam.o get_contam_c2f.o js_kdtree.o -o get_contam.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c get_ptcl.f90 -o get_ptcl.o
gcc -fopenmp -fPIC -mcmodel=large -c get_ptcl_c2f.c -o get_ptcl_c2f.o
gcc -shared -fopenmp get_ptcl.o get_ptcl_c2f.o -o get_ptcl.so -lgfortran

gfortran -fopenmp -fPIC -c prop_time.f90 -o prop_time.o
gcc -fopenmp -fPIC -c prop_time_c2f.c -o prop_time_c2f.o
gcc -shared -fopenmp prop_time.o prop_time_c2f.o -o prop_time.so -lgfortran

gfortran -fopenmp -fPIC -c jsamr2cell_totnum.f90 -o jsamr2cell_totnum.o
gcc -fopenmp -fPIC -c jsamr2cell_totnum_c2f.c -o jsamr2cell_totnum_c2f.o
gcc -shared -fopenmp jsamr2cell_totnum.o jsamr2cell_totnum_c2f.o -o jsamr2cell_totnum.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c jsamr2cell.f90 -o jsamr2cell.o
gcc -fopenmp -fPIC -mcmodel=large -c jsamr2cell_c2f.c -o jsamr2cell_c2f.o
gcc -shared -fopenmp -mcmodel=large jsamr2cell.o jsamr2cell_c2f.o -o jsamr2cell.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c js_getpt_ft.f90 -o js_getpt_ft.o
gcc -fopenmp -fPIC -mcmodel=large -c js_getpt_ft_c2f.c -o js_getpt_ft_c2f.o
gcc -shared -fopenmp -mcmodel=large js_getpt_ft.o js_getpt_ft_c2f.o js_kdtree.o -o js_getpt_ft.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -fno-range-check -c rv_match.f90 -o rv_match.o
gcc -fopenmp -fPIC -mcmodel=large -c rv_match_c2f.c -o rv_match_c2f.o
gcc -shared -fopenmp rv_match.o rv_match_c2f.o -o rv_match.so -lgfortran

gfortran -fopenmp -fPIC -c prop_sfr.f90 -o prop_sfr.o
gcc -fopenmp -fPIC -c prop_sfr_c2f.c -o prop_sfr_c2f.o
gcc -shared -fopenmp prop_sfr.o prop_sfr_c2f.o -o prop_sfr.so
chmod 777 *.o *.so

gfortran -fopenmp -fPIC -mcmodel=large -c js_kde_2dhisto.f90 -o js_kde_2dhisto.o
gcc -fopenmp -fPIC -mcmodel=large -c js_kde_2dhisto_c2f.c -o js_kde_2dhisto_c2f.o
gcc -shared -fopenmp -mcmodel=large js_kde_2dhisto.o js_kde_2dhisto_c2f.o -o js_kde_2dhisto.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c js_kde_gauss_binning.f90 -o js_kde_gauss_binning.o
gcc -fopenmp -fPIC -mcmodel=large -c js_kde_gauss_binning_c2f.c -o js_kde_gauss_binning_c2f.o
gcc -shared -fopenmp -mcmodel=large js_kde_gauss_binning.o js_kde_gauss_binning_c2f.o -o js_kde_gauss_binning.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c js_kde_gauss_pts.f90 -o js_kde_gauss_pts.o
gcc -fopenmp -fPIC -mcmodel=large -c js_kde_gauss_pts_c2f.c -o js_kde_gauss_pts_c2f.o
gcc -shared -fopenmp -mcmodel=large js_kde_gauss_pts.o js_kde_gauss_pts_c2f.o -o js_kde_gauss_pts.so -lgfortran

gfortran -fopenmp -fPIC -mcmodel=large -c js_gasmap.f90 -o js_gasmap.o
gcc -fopenmp -fPIC -mcmodel=large -c js_gasmap_c2f.c -o js_gasmap_c2f.o
gcc -shared -fopenmp -mcmodel=large js_gasmap.o js_gasmap_c2f.o fft/fftlib.o -o js_gasmap.so -lgfortran

