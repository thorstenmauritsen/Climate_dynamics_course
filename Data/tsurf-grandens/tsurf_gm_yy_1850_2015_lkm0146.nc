CDF   �   
      lon       lat       time       bnds            CDI       ?Climate Data Interface version 1.6.9 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      Fri Oct 02 15:31:38 2015: cdo -s -yearmonmean -fldsum -mul -gridweights /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/echam6/lkm0146_tsurf_gm_yy_1850_2015.nc
Fri Oct 02 15:31:35 2015: cdo -s -mergetime /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2005.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_2006_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc
Fri Oct 02 15:31:29 2015: cdo -f nc -t echam6 -s -r -select,name=tsurf /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2006.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2007.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2008.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2009.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2010.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2011.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2012.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45146/outdata/echam6/jkrcp45146_echam6_BOT_mm_2013.grb      source        ECHAM6     institution       $Max-Planck-Institute for Meteorology   CDO       @Climate Data Operators version 1.6.9 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X           
�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           
�   time               standard_name         time   bounds        	time_bnds      units         hours since 1850-01-31 23:52:00    calendar      proleptic_gregorian         
�   	time_bnds                            
�   tsurf                         	long_name         surface temperature    units         K      code         �   table            �   
_FillValue        ����   missing_value         ����        
�                @�             @�P     C���@��     @�     @�b     C��@@�     @�     @��     C��<@ݑ     @ٰ     @��     C��@��    @�     @�	     C��	@�V�    @�f     @�P     C��_@렀    @�     @�     C���@��    @��     @���    C��2@�@    @�     @�     C��@�:�    @�B�    @�7�    C��7@�_�    @�f     @�\�    C���@��@    @��     @��     C��(@���    @���    @���    C���@��@    @��     @��     C�ە@��@    @���    @��     C�ϲA �`    A @    A�    C���A�     A     A�    C��:A��    A0�    A+@    C��KA�`    AB�    A=�    C���A�     AU     AO�    C��A��    Af�    Aa@    C��TA��    Ax�    As     C��~A     A�@    A��    C��A	�    A��    A	�@    C���A
*�    A	��    A
�     C��|A<`    A
�@    A��    C��AN�    A�     A�@    C��A`�    A�    A�     C�ٯAr`    A�@    A��    C��A�     A     A@    C��JAKP    A�    A��    C��$A�0    A�     A`    C���A]    A     A�@    C��2A��    A��    A%     C��EAo0    A0�    A�`    C���A�    A�     A7@    C��lA��    AB�    A�     C��mA	�    A��    AI     C��A�    AT�    A�@    C��{A�    A��    A[     C��GA��    Af�    A�     C���A-�    A�    Al�    C���A��    Ax�    A�     C��gA?�    A�    A     C��NAȰ    A��    A�    C���AQ�    A�    A��    C��A��    A�`    A     C��VAc�    A%�    A��    C��SA�    A��    A+�    C��IAup    A7`    A��    C��ZA�P    A�@    A=�    C���A�0    AI     A�`    C��A    A�     AO@    C���A��    AZ�    A�     C���A"0    A��    Aa`    C�İA�    Am     A�@    C���A3�    A��    As     C��A��    A~�    A�     C���AF    A�    A�@    C��%A��    A��    A     C���A +�    A �    A K�    C���A pX    A QP    A ��    C��!A ��    A ��    A Ԑ    C��HA �h    A �`    A!     C���A!=�    A!�    A!]p    C�ةA!�H    A!c@    A!��    C���A!��    A!��    A!�    C���A"X    A!�P    A"*�    C��A"O�    A"0�    A"o`    C��aA"�8    A"u0    A"��    C���A"��    A"��    A"�p    C�ΙA#H    A"�@    A#<�    C��A#a�    A#B�    A#�P    C��A#�(    A#�     A#��    C��A#��    A#ː    A$
`    C���A$/8    A$0    A$N�    C�� A$s�    A$T�    A$�@    C���A$�    A$�    A$װ    C��A$��    A$݀    A%P    C���A%A(    A%"     A%`�    C���A%��    A%f�    A%�0    C��VA%�    A%�     A%�    C��A&�    A%�p    A&.@    C���A&S    A&4    A&r�    C��iA&��    A&x�    A&�     C��zA&��    A&��    A&��    C�ܠA' �    A'`    A'@0    C��A'e    A'F     A'��    C��jA'�x    A'�p    A'�    C��-A'��    A'��    A(�    C��A(2�    A(P    A(R     C�A(v�    A(W�    A(��    C�
�A(�h    A(�`    A(�     C���A(��    A(��    A)p    C���A)Dx    A)%@    A)d    C���A)��    A)i�    A)��    C��A)�X    A)�P    A)��    C��^A*�    A)��    A*1`    C��:A*Vh    A*70    A*v     C���A*��    A*{�    A*�p    C���A*�H    A*�@    A*��    C��,A+#�    A+�    A+CP    C��A+hX    A+I     A+��    C���A+��    A+��    A+�`    C���A+�8    A+�0    A,�    C��A,5�    A,�    A,U@    C���A,zH    A,[    A,��    C���A,��    A,��    A,�P    C�ܷA-(    A,�     A-"�    C��A-G�    A-(�    A-g0    C���A-�8    A-m     A-��    C�A-Ш    A-��    A-�@    C���A.    A-�    A.4�    C��{A.Y�    A.:�    A.y     C��YA.�(    A.~�    A.��    C���A.�    A.Ð    A/0    C��A/'    A/     A/F�    C��A/kx    A/Lp    A/�    C���A/�    A/��    A/ϰ    C���A/�    A/Հ    A0
    C��iA0|    A0�    A0,H    C��{A0>�    A0/0    A0N�    C��}A0a    A0Qh    A0p�    C��A0�<    A0s�    A0�    C���A0�t    A0��    A0�@    C��A0Ǭ    A0�(    A0�x    C��EA0��    A0�`    A0��    C���A14    A0��    A1     C��A1.l    A1�    A1>8    C��=A1P�    A1A     A1`p    C��HA1r�    A1cX    A1��    C�WA1�,    A1��    A1��    C�"A1�d    A1��    A1�0    C�eA1ٜ    A1�    A1�h    C���A1��    A1�P    A2�    C��>A2$    A2�    A2-�    C��A2@\    A20�    A2P(    C�=A2b�    A2S    A2r`    C��A2��    A2uH    A2��    C�%PA2�    A2��    A2��    C�"A2�T    A2��    A2�     C�"A2�    A2�    A2�X    C��A3�    A2�@    A3�    C��AA30    A3 �    A3?�    C���A3RL    A3B�    A3b    C��A3t�    A3e     A3�P    C��A3��    A3�8    A3��    C��A3�    A3��    A3��    C�%A3�D    A3��    A3�    C�!A3�|    A3��    A4H    C�/�A4�    A40    A4/�    C�.�A4B    A42�    A4Q�    C�<�A4d<    A4T�    A4t    C�2�A4�t    A4v�    A4�@    C�05A4��    A4�(    A4��    C�<1A4��    A4�x    A4��    C�>A4�4    A4ݰ    A4�     C�> A5l    A4��    A58    C�I�A51�    A5"     A5A�    C�GA5S�    A5Dp    A5c�    C�D�A5v,    A5f�    A5��    C�B�A5�d    A5��    A5�0    C�F�A5��    A5�    A5ʀ    C�N�A5��    A5�h    A5�    C�W�A5�$    A5�    A6�    C�W�A6!\    A6�    A61(    C�J�