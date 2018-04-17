CDF   �   
      lon       lat       time       bnds            CDI       ?Climate Data Interface version 1.6.9 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      Fri Oct 02 15:35:44 2015: cdo -s -yearmonmean -fldsum -mul -gridweights /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/echam6/lkm0148_tsurf_gm_yy_1850_2015.nc
Fri Oct 02 15:35:40 2015: cdo -s -mergetime /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2005.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_2006_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc
Fri Oct 02 15:35:34 2015: cdo -f nc -t echam6 -s -r -select,name=tsurf /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2006.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2007.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2008.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2009.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2010.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2011.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2012.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45148/outdata/echam6/jkrcp45148_echam6_BOT_mm_2013.grb      source        ECHAM6     institution       $Max-Planck-Institute for Meteorology   CDO       @Climate Data Operators version 1.6.9 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X           
�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           
�   time               standard_name         time   bounds        	time_bnds      units         hours since 1850-01-31 23:52:00    calendar      proleptic_gregorian         
�   	time_bnds                            
�   tsurf                         	long_name         surface temperature    units         K      code         �   table            �   
_FillValue        ����   missing_value         ����        
�                @�             @�P     C��a@��     @�     @�b     C��@�     @�     @��     C�܉@ݑ     @ٰ     @��     C���@��    @�     @�	     C��*@�V�    @�f     @�P     C���@렀    @�     @�     C���@��    @��     @���    C��@�@    @�     @�     C���@�:�    @�B�    @�7�    C���@�_�    @�f     @�\�    C��`@��@    @��     @��     C�ݻ@���    @���    @���    C���@��@    @��     @��     C�̍@��@    @���    @��     C��<A �`    A @    A�    C��oA�     A     A�    C��oA��    A0�    A+@    C��bA�`    AB�    A=�    C���A�     AU     AO�    C���A��    Af�    Aa@    C���A��    Ax�    As     C��xA     A�@    A��    C��BA	�    A��    A	�@    C���A
*�    A	��    A
�     C��<A<`    A
�@    A��    C���AN�    A�     A�@    C���A`�    A�    A�     C��~Ar`    A�@    A��    C��?A�     A     A@    C��AKP    A�    A��    C�� A�0    A�     A`    C��:A]    A     A�@    C���A��    A��    A%     C��hAo0    A0�    A�`    C��A�    A�     A7@    C�ѭA��    AB�    A�     C���A	�    A��    AI     C��cA�    AT�    A�@    C��BA�    A��    A[     C���A��    Af�    A�     C���A-�    A�    Al�    C���A��    Ax�    A�     C���A?�    A�    A     C���AȰ    A��    A�    C��.AQ�    A�    A��    C�ãA��    A�`    A     C�ċAc�    A%�    A��    C�ˆA�    A��    A+�    C��)Aup    A7`    A��    C���A�P    A�@    A=�    C��gA�0    AI     A�`    C�БA    A�     AO@    C�̃A��    AZ�    A�     C��cA"0    A��    Aa`    C���A�    Am     A�@    C���A3�    A��    As     C���A��    A~�    A�     C��YAF    A�    A�@    C��TA��    A��    A     C���A +�    A �    A K�    C��bA pX    A QP    A ��    C�ڱA ��    A ��    A Ԑ    C��A �h    A �`    A!     C��A!=�    A!�    A!]p    C���A!�H    A!c@    A!��    C��A!��    A!��    A!�    C��zA"X    A!�P    A"*�    C��dA"O�    A"0�    A"o`    C��kA"�8    A"u0    A"��    C��}A"��    A"��    A"�p    C��A#H    A"�@    A#<�    C��2A#a�    A#B�    A#�P    C���A#�(    A#�     A#��    C��A#��    A#ː    A$
`    C��{A$/8    A$0    A$N�    C��A$s�    A$T�    A$�@    C��A$�    A$�    A$װ    C�iA$��    A$݀    A%P    C��tA%A(    A%"     A%`�    C���A%��    A%f�    A%�0    C���A%�    A%�     A%�    C�ӪA&�    A%�p    A&.@    C�ֱA&S    A&4    A&r�    C�ܲA&��    A&x�    A&�     C���A&��    A&��    A&��    C�̧A' �    A'`    A'@0    C��%A'e    A'F     A'��    C��`A'�x    A'�p    A'�    C��6A'��    A'��    A(�    C��aA(2�    A(P    A(R     C��cA(v�    A(W�    A(��    C��NA(�h    A(�`    A(�     C�׳A(��    A(��    A)p    C��A)Dx    A)%@    A)d    C��A)��    A)i�    A)��    C��wA)�X    A)�P    A)��    C��PA*�    A)��    A*1`    C���A*Vh    A*70    A*v     C�ޣA*��    A*{�    A*�p    C��WA*�H    A*�@    A*��    C���A+#�    A+�    A+CP    C��A+hX    A+I     A+��    C��A+��    A+��    A+�`    C��|A+�8    A+�0    A,�    C�ܤA,5�    A,�    A,U@    C��DA,zH    A,[    A,��    C���A,��    A,��    A,�P    C���A-(    A,�     A-"�    C��A-G�    A-(�    A-g0    C���A-�8    A-m     A-��    C��ZA-Ш    A-��    A-�@    C��qA.    A-�    A.4�    C���A.Y�    A.:�    A.y     C�#�A.�(    A.~�    A.��    C��A.�    A.Ð    A/0    C��oA/'    A/     A/F�    C��MA/kx    A/Lp    A/�    C���A/�    A/��    A/ϰ    C��5A/�    A/Հ    A0
    C���A0|    A0�    A0,H    C��A0>�    A0/0    A0N�    C�
+A0a    A0Qh    A0p�    C��A0�<    A0s�    A0�    C��A0�t    A0��    A0�@    C��A0Ǭ    A0�(    A0�x    C�A0��    A0�`    A0��    C��1A14    A0��    A1     C��MA1.l    A1�    A1>8    C�9A1P�    A1A     A1`p    C�JA1r�    A1cX    A1��    C��A1�,    A1��    A1��    C�!MA1�d    A1��    A1�0    C�(�A1ٜ    A1�    A1�h    C�FA1��    A1�P    A2�    C�HA2$    A2�    A2-�    C�2A2@\    A20�    A2P(    C�UA2b�    A2S    A2r`    C�|A2��    A2uH    A2��    C�yA2�    A2��    A2��    C�/A2�T    A2��    A2�     C�A2�    A2�    A2�X    C�!^A3�    A2�@    A3�    C���A30    A3 �    A3?�    C���A3RL    A3B�    A3b    C�hA3t�    A3e     A3�P    C�FA3��    A3�8    A3��    C�,�A3�    A3��    A3��    C�*�A3�D    A3��    A3�    C��A3�|    A3��    A4H    C�"sA4�    A40    A4/�    C�$�A4B    A42�    A4Q�    C�5�A4d<    A4T�    A4t    C�3hA4�t    A4v�    A4�@    C�G�A4��    A4�(    A4��    C�]-A4��    A4�x    A4��    C�KA4�4    A4ݰ    A4�     C�9YA5l    A4��    A58    C�-�A51�    A5"     A5A�    C��A5S�    A5Dp    A5c�    C�(vA5v,    A5f�    A5��    C�O	A5�d    A5��    A5�0    C�D�A5��    A5�    A5ʀ    C�;�A5��    A5�h    A5�    C�<�A5�$    A5�    A6�    C�YA6!\    A6�    A61(    C�Z�