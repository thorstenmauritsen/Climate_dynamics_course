CDF   �   
      lon       lat       time       bnds            CDI       ?Climate Data Interface version 1.6.9 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      Fri Oct 02 14:04:56 2015: cdo -s -yearmonmean -fldsum -mul -gridweights /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/echam6/lkm0108_tsurf_gm_yy_1850_2015.nc
Fri Oct 02 14:04:53 2015: cdo -s -mergetime /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2005.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_2006_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc
Fri Oct 02 14:04:46 2015: cdo -f nc -t echam6 -s -r -select,name=tsurf /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2006.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2007.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2008.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2009.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2010.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2011.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2012.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45108/outdata/echam6/jkrcp45108_echam6_BOT_mm_2013.grb      source        ECHAM6     institution       $Max-Planck-Institute for Meteorology   CDO       @Climate Data Operators version 1.6.9 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X           
�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           
�   time               standard_name         time   bounds        	time_bnds      units         hours since 1850-01-31 23:52:00    calendar      proleptic_gregorian         
�   	time_bnds                            
�   tsurf                         	long_name         surface temperature    units         K      code         �   table            �   
_FillValue        ����   missing_value         ����        
�                @�             @�P     C�͓@��     @�     @�b     C�� @�     @�     @��     C��?@ݑ     @ٰ     @��     C���@��    @�     @�	     C�˅@�V�    @�f     @�P     C���@렀    @�     @�     C��7@��    @��     @���    C��@�@    @�     @�     C��7@�:�    @�B�    @�7�    C��t@�_�    @�f     @�\�    C���@��@    @��     @��     C��@���    @���    @���    C���@��@    @��     @��     C��b@��@    @���    @��     C���A �`    A @    A�    C�ʒA�     A     A�    C��A��    A0�    A+@    C��5A�`    AB�    A=�    C��'A�     AU     AO�    C���A��    Af�    Aa@    C���A��    Ax�    As     C��
A     A�@    A��    C���A	�    A��    A	�@    C��A
*�    A	��    A
�     C��A<`    A
�@    A��    C��dAN�    A�     A�@    C��rA`�    A�    A�     C��JAr`    A�@    A��    C��/A�     A     A@    C���AKP    A�    A��    C���A�0    A�     A`    C��;A]    A     A�@    C��"A��    A��    A%     C��Ao0    A0�    A�`    C��DA�    A�     A7@    C���A��    AB�    A�     C��,A	�    A��    AI     C���A�    AT�    A�@    C��A�    A��    A[     C���A��    Af�    A�     C��A-�    A�    Al�    C���A��    Ax�    A�     C��vA?�    A�    A     C��AȰ    A��    A�    C��OAQ�    A�    A��    C��tA��    A�`    A     C��KAc�    A%�    A��    C���A�    A��    A+�    C��GAup    A7`    A��    C���A�P    A�@    A=�    C���A�0    AI     A�`    C���A    A�     AO@    C�˰A��    AZ�    A�     C��,A"0    A��    Aa`    C��kA�    Am     A�@    C���A3�    A��    As     C��A��    A~�    A�     C��fAF    A�    A�@    C��:A��    A��    A     C���A +�    A �    A K�    C��wA pX    A QP    A ��    C��=A ��    A ��    A Ԑ    C���A �h    A �`    A!     C���A!=�    A!�    A!]p    C��!A!�H    A!c@    A!��    C��AA!��    A!��    A!�    C��A"X    A!�P    A"*�    C��	A"O�    A"0�    A"o`    C�ͷA"�8    A"u0    A"��    C��PA"��    A"��    A"�p    C��A#H    A"�@    A#<�    C���A#a�    A#B�    A#�P    C��IA#�(    A#�     A#��    C��+A#��    A#ː    A$
`    C�ɬA$/8    A$0    A$N�    C��[A$s�    A$T�    A$�@    C���A$�    A$�    A$װ    C��cA$��    A$݀    A%P    C��'A%A(    A%"     A%`�    C���A%��    A%f�    A%�0    C��A%�    A%�     A%�    C��A&�    A%�p    A&.@    C�ƭA&S    A&4    A&r�    C��LA&��    A&x�    A&�     C��	A&��    A&��    A&��    C��hA' �    A'`    A'@0    C��cA'e    A'F     A'��    C�̀A'�x    A'�p    A'�    C���A'��    A'��    A(�    C��{A(2�    A(P    A(R     C���A(v�    A(W�    A(��    C���A(�h    A(�`    A(�     C��7A(��    A(��    A)p    C��TA)Dx    A)%@    A)d    C��KA)��    A)i�    A)��    C��QA)�X    A)�P    A)��    C���A*�    A)��    A*1`    C��7A*Vh    A*70    A*v     C��A*��    A*{�    A*�p    C��mA*�H    A*�@    A*��    C��uA+#�    A+�    A+CP    C��A+hX    A+I     A+��    C��CA+��    A+��    A+�`    C���A+�8    A+�0    A,�    C��IA,5�    A,�    A,U@    C��OA,zH    A,[    A,��    C���A,��    A,��    A,�P    C���A-(    A,�     A-"�    C���A-G�    A-(�    A-g0    C��A-�8    A-m     A-��    C�3A-Ш    A-��    A-�@    C��A.    A-�    A.4�    C��A.Y�    A.:�    A.y     C��uA.�(    A.~�    A.��    C���A.�    A.Ð    A/0    C�ͰA/'    A/     A/F�    C��`A/kx    A/Lp    A/�    C���A/�    A/��    A/ϰ    C��'A/�    A/Հ    A0
    C���A0|    A0�    A0,H    C��iA0>�    A0/0    A0N�    C��CA0a    A0Qh    A0p�    C��"A0�<    A0s�    A0�    C���A0�t    A0��    A0�@    C��`A0Ǭ    A0�(    A0�x    C�A0��    A0�`    A0��    C�A14    A0��    A1     C��|A1.l    A1�    A1>8    C���A1P�    A1A     A1`p    C�jA1r�    A1cX    A1��    C��A1�,    A1��    A1��    C��A1�d    A1��    A1�0    C��A1ٜ    A1�    A1�h    C��@A1��    A1�P    A2�    C�A2$    A2�    A2-�    C�+#A2@\    A20�    A2P(    C� 4A2b�    A2S    A2r`    C��A2��    A2uH    A2��    C��A2�    A2��    A2��    C�aA2�T    A2��    A2�     C�
�A2�    A2�    A2�X    C��A3�    A2�@    A3�    C�[A30    A3 �    A3?�    C��A3RL    A3B�    A3b    C�	A3t�    A3e     A3�P    C���A3��    A3�8    A3��    C��A3�    A3��    A3��    C�A3�D    A3��    A3�    C��A3�|    A3��    A4H    C�'+A4�    A40    A4/�    C�0�A4B    A42�    A4Q�    C�A�A4d<    A4T�    A4t    C�L�A4�t    A4v�    A4�@    C�<�A4��    A4�(    A4��    C�5A4��    A4�x    A4��    C�D
A4�4    A4ݰ    A4�     C�=�A5l    A4��    A58    C�@HA51�    A5"     A5A�    C�2�A5S�    A5Dp    A5c�    C�8�A5v,    A5f�    A5��    C�:�A5�d    A5��    A5�0    C�?nA5��    A5�    A5ʀ    C�qUA5��    A5�h    A5�    C�_;A5�$    A5�    A6�    C�J�A6!\    A6�    A61(    C�I�