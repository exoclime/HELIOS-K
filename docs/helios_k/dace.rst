Dace Database
=============


The `DACE Opacity Database` (https://dace.unige.ch) contain a large number of calculated Opacities. Most species from the EXOMOL, HITRAN, NIST, VALD, and Kurucz databases
are available there. 

This section shows how to calculate new opacities in the same format as in DACE.


Step 1: Download and prepare the linelists.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the instructions on the section about the different linelist databse :ref:`exomol` :ref:`hitran`, :ref:`hitemp`, :ref:`kurucz`, :ref:`nist`, and :ref:`vald`.
This step will create the binary line lists and the species.param files needed to run Helios-k 



Step 2: Create a script to run Helios-k for all desired T P and nu values.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The values for the tmeperature, pressure and wavenumbers in DACE are set in the following way:

* wavenumber :math:`\nu`: from zero to the maximum available wavenuumber of the individual line list. Round up the maximum to the next step of 1000 :math:`cm^{-1}`.
* Temperature T: Check the maximum available temperature in the corresponding partition function file. The DACE database uses the following temperatures::


	50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, (step = 50K)
	700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500,	   (step = 100K)
	1700, 1900, 2100, 2300, 2500, 2700, 2900, ...              (step = 200K)


* Pressure p: The DACE database uses the following pressure values wiht the corresponding notation::

	n800	:math:`10^{-8.00} bar = 9.8692326700000007768e-09 atm
	n766	:math:`10^{-7.66} bar = 2.1262617228243976416e-08 atm
	n733	:math:`10^{-7.33} bar = 4.5808920157398465827e-08 atm
	n700	:math:`10^{-7.00} bar = 9.8692326699999991224e-08 atm
	n666	:math:`10^{-6.66} bar = 2.1262617228243975093e-07 atm
	n633	:math:`10^{-6.33} bar = 4.580892015739846318e-07 atm
	n600	:math:`10^{-6.00} bar = 9.8692326700000004459e-07 atm
	n566	:math:`10^{-5.66} bar = 2.1262617228243974034e-06 atm
	n533	:math:`10^{-5.33} bar = 4.5808920157398460004e-06 atm
	n500	:math:`10^{-5.00} bar = 9.8692326700000012929e-06 atm
	n466	:math:`10^{-4.66} bar = 2.1262617228243976575e-05 atm
	n433	:math:`10^{-4.33} bar = 4.5808920157398456616e-05 atm
	n400	:math:`10^{-4.00} bar = 9.8692326700000002765e-05 atm
	n366	:math:`10^{-3.66} bar = 0.00021262617228243998937 atm
	n333	:math:`10^{-3.33} bar = 0.00045808920157398413248 atm
	n300	:math:`10^{-3.00} bar = 0.00098692326700000002765 atm
	n266	:math:`10^{-2.66} bar = 0.0021262617228243996768 atm
	n233	:math:`10^{-2.33} bar = 0.0045808920157398413248 atm
	n200	:math:`10^{-2.00} bar = 0.0098692326700000002765 atm
	n166	:math:`10^{-1.66} bar = 0.021262617228243983758 atm
	n133	:math:`10^{-1.33} bar = 0.045808920157398416717 atm
	n100	:math:`10^{-1.00} bar = 0.098692326699999999295 atm
	n066	:math:`10^{-0.66} bar = 0.21262617228243987921 atm
	n033	:math:`10^{-0.33} bar = 0.45808920157398430595 atm
	p000	:math:`10^{0.00} bar = 0.98692326699999999295 atm
	p033	:math:`10^{0.33} bar = 2.1262617228243989587 atm
	p066	:math:`10^{0.66} bar = 4.5808920157398427264 atm
	p100	:math:`10^{1.00} bar = 9.8692326700000005957 atm
	p133	:math:`10^{1.33} bar = 21.262617228243986034 atm
	p166	:math:`10^{1.66} bar = 45.808920157398411277 atm
	p200	:math:`10^{2.00} bar = 98.692326699999995299 atm
	p233	:math:`10^{2.33} bar = 212.62617228243999534 atm
	p266	:math:`10^{2.66} bar = 458.0892015739841554 atm
	p300	:math:`10^{3.00} bar = 986.92326700000000983 atm

	




This can be done with the initialPressure code (availalble in the repository tools directory).
Choose the maximum temperature according to the maximum value in the partition function file::

        ./initialPresure -M 27_1_hit -nMax 3 -TMax 38


Step 3: Create a new folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a new folder  with the name of the species e.g.  ``12C2-1H6__HITRAN_e2``


The ending indicates a resolution of 10^-2 in binary format


Step 4: Copy files
~~~~~~~~~~~~~~~~~~

Copy the run sript, heliosk, and param.dat into that folder


Run script::

        bash script<name>

This will run Helios-K for all necessary input parameters


Step 5:
~~~~~~~

        mkdir 12C2-1H6__HITRAN_e0Dat
        mkdir 12C2-1H6__HITRAN_e2b

        mkdir 12C2-1H6__HITRAN_e0Dat

        ./reduce -b 0 -reduce 100 -dnu 0.01 -numin 0 -numax 4000 -Tmin 0 -Tmax 4900 -indir 12C2-1H6__HITRAN_e2 -outdir 12C2-1H6__HITRAN_e0Dat

        ./reduce -b 1 -reduce 1 -dnu 0.01 -numin 0 -numax 4000 -Tmin 0 -Tmax 4900 -indir 12C2-1H6__HITRAN_e2 -outdir 12C2-1H6__HITRAN_e2b

        ./reduce -b 0 -reduce 100 -dnu 0.01 -numin 0 -numax 1000 -Tmin 0 -Tmax 4900 -indir 12C-1H3-13C-1H3__HITRAN_e2 -outdir 12C-1H3-13C-1H3__HITRAN_e0Dat
        ./reduce -b 1 -reduce 1 -dnu 0.01 -numin 0 -numax 1000 -Tmin 0 -Tmax 4900 -indir 12C-1H3-13C-1H3__HITRANe2 -outdir 12C-1H3-13C-1H3__HITRAN_e2b


