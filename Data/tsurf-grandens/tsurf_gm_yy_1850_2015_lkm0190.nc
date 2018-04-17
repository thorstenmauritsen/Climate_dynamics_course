CDF   �   
      lon       lat       time       bnds            CDI       ?Climate Data Interface version 1.6.9 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      Fri Oct 02 17:08:14 2015: cdo -s -yearmonmean -fldsum -mul -gridweights /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/echam6/lkm0190_tsurf_gm_yy_1850_2015.nc
Fri Oct 02 17:08:12 2015: cdo -s -mergetime /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2005.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_2006_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc
Fri Oct 02 17:08:04 2015: cdo -f nc -t echam6 -s -r -select,name=tsurf /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2006.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2007.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2008.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2009.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2010.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2011.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2012.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45190/outdata/echam6/jkrcp45190_echam6_BOT_mm_2013.grb      source        ECHAM6     institution       $Max-Planck-Institute for Meteorology   CDO       @Climate Data Operators version 1.6.9 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X           
�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           
�   time               standard_name         time   bounds        	time_bnds      units         hours since 1850-01-31 23:52:00    calendar      proleptic_gregorian         
�   	time_bnds                            
�   tsurf                         	long_name         surface temperature    units         K      code         �   table            �   
_FillValue        ����   missing_value         ����        
�                @�             @�P     C��n@��     @�     @�b     C��@�     @�     @��     C���@ݑ     @ٰ     @��     C��@��    @�     @�	     C���@�V�    @�f     @�P     C�� @렀    @�     @�     C��,@��    @��     @���    C���@�@    @�     @�     C��$@�:�    @�B�    @�7�    C��L@�_�    @�f     @�\�    C��@��@    @��     @��     C��1@���    @���    @���    C���@��@    @��     @��     C���@��@    @���    @��     C��jA �`    A @    A�    C��1A�     A     A�    C���A��    A0�    A+@    C���A�`    AB�    A=�    C�A�     AU     AO�    C��CA��    Af�    Aa@    C��;A��    Ax�    As     C��A     A�@    A��    C��A	�    A��    A	�@    C��BA
*�    A	��    A
�     C��KA<`    A
�@    A��    C��/AN�    A�     A�@    C��eA`�    A�    A�     C��Ar`    A�@    A��    C��gA�     A     A@    C���AKP    A�    A��    C��A�0    A�     A`    C��A]    A     A�@    C��`A��    A��    A%     C��Ao0    A0�    A�`    C���A�    A�     A7@    C���A��    AB�    A�     C�ƎA	�    A��    AI     C��
A�    AT�    A�@    C��<A�    A��    A[     C�� A��    Af�    A�     C��mA-�    A�    Al�    C��lA��    Ax�    A�     C���A?�    A�    A     C���AȰ    A��    A�    C���AQ�    A�    A��    C��2A��    A�`    A     C��oAc�    A%�    A��    C��nA�    A��    A+�    C���Aup    A7`    A��    C��|A�P    A�@    A=�    C��A�0    AI     A�`    C��A    A�     AO@    C��nA��    AZ�    A�     C���A"0    A��    Aa`    C���A�    Am     A�@    C���A3�    A��    As     C���A��    A~�    A�     C��zAF    A�    A�@    C��A��    A��    A     C��MA +�    A �    A K�    C��#A pX    A QP    A ��    C���A ��    A ��    A Ԑ    C��SA �h    A �`    A!     C�̶A!=�    A!�    A!]p    C���A!�H    A!c@    A!��    C���A!��    A!��    A!�    C�ԨA"X    A!�P    A"*�    C��A"O�    A"0�    A"o`    C��-A"�8    A"u0    A"��    C��YA"��    A"��    A"�p    C��?A#H    A"�@    A#<�    C�ǺA#a�    A#B�    A#�P    C�ԬA#�(    A#�     A#��    C��rA#��    A#ː    A$
`    C���A$/8    A$0    A$N�    C���A$s�    A$T�    A$�@    C��A$�    A$�    A$װ    C��<A$��    A$݀    A%P    C���A%A(    A%"     A%`�    C��,A%��    A%f�    A%�0    C���A%�    A%�     A%�    C��	A&�    A%�p    A&.@    C���A&S    A&4    A&r�    C���A&��    A&x�    A&�     C��<A&��    A&��    A&��    C��5A' �    A'`    A'@0    C���A'e    A'F     A'��    C��A'�x    A'�p    A'�    C��0A'��    A'��    A(�    C��KA(2�    A(P    A(R     C���A(v�    A(W�    A(��    C�
cA(�h    A(�`    A(�     C��6A(��    A(��    A)p    C��A)Dx    A)%@    A)d    C���A)��    A)i�    A)��    C��CA)�X    A)�P    A)��    C��xA*�    A)��    A*1`    C��A*Vh    A*70    A*v     C���A*��    A*{�    A*�p    C���A*�H    A*�@    A*��    C��A+#�    A+�    A+CP    C���A+hX    A+I     A+��    C��A+��    A+��    A+�`    C���A+�8    A+�0    A,�    C��vA,5�    A,�    A,U@    C���A,zH    A,[    A,��    C��A,��    A,��    A,�P    C��A-(    A,�     A-"�    C��6A-G�    A-(�    A-g0    C���A-�8    A-m     A-��    C���A-Ш    A-��    A-�@    C��A.    A-�    A.4�    C��HA.Y�    A.:�    A.y     C��A.�(    A.~�    A.��    C���A.�    A.Ð    A/0    C���A/'    A/     A/F�    C���A/kx    A/Lp    A/�    C���A/�    A/��    A/ϰ    C��A/�    A/Հ    A0
    C��A0|    A0�    A0,H    C��A0>�    A0/0    A0N�    C��A0a    A0Qh    A0p�    C���A0�<    A0s�    A0�    C��A0�t    A0��    A0�@    C�6RA0Ǭ    A0�(    A0�x    C�+�A0��    A0�`    A0��    C��DA14    A0��    A1     C���A1.l    A1�    A1>8    C���A1P�    A1A     A1`p    C���A1r�    A1cX    A1��    C�rA1�,    A1��    A1��    C�]A1�d    A1��    A1�0    C��A1ٜ    A1�    A1�h    C��A1��    A1�P    A2�    C�-A2$    A2�    A2-�    C���A2@\    A20�    A2P(    C���A2b�    A2S    A2r`    C�	A2��    A2uH    A2��    C��A2�    A2��    A2��    C�-�A2�T    A2��    A2�     C�C�A2�    A2�    A2�X    C�"�A3�    A2�@    A3�    C�ۏA30    A3 �    A3?�    C���A3RL    A3B�    A3b    C�)A3t�    A3e     A3�P    C�A3��    A3�8    A3��    C��A3�    A3��    A3��    C��A3�D    A3��    A3�    C�@�A3�|    A3��    A4H    C�CA4�    A40    A4/�    C�&eA4B    A42�    A4Q�    C�0�A4d<    A4T�    A4t    C�:�A4�t    A4v�    A4�@    C�K(A4��    A4�(    A4��    C�<lA4��    A4�x    A4��    C�7A4�4    A4ݰ    A4�     C��A5l    A4��    A58    C�A51�    A5"     A5A�    C�9cA5S�    A5Dp    A5c�    C�I�A5v,    A5f�    A5��    C�^;A5�d    A5��    A5�0    C�J�A5��    A5�    A5ʀ    C�8RA5��    A5�h    A5�    C�H�A5�$    A5�    A6�    C�[�A6!\    A6�    A61(    C�F�