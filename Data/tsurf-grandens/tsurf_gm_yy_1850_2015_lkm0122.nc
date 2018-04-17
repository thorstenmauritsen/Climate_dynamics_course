CDF   �   
      lon       lat       time       bnds            CDI       ?Climate Data Interface version 1.6.9 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      Fri Oct 02 14:37:27 2015: cdo -s -yearmonmean -fldsum -mul -gridweights /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/echam6/lkm0122_tsurf_gm_yy_1850_2015.nc
Fri Oct 02 14:37:24 2015: cdo -s -mergetime /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2005.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_2006_2015.nc /scratch/mpi/mpioes/m300297/paper-01/raw-data/tmp/tmp_tsurf_1850_2015.nc
Fri Oct 02 14:37:18 2015: cdo -f nc -t echam6 -s -r -select,name=tsurf /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2006.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2007.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2008.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2009.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2010.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2011.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2012.grb /scratch/mpi/mpioes/mh0033/m300297/hist100_2006-2015/jkrcp45122/outdata/echam6/jkrcp45122_echam6_BOT_mm_2013.grb      source        ECHAM6     institution       $Max-Planck-Institute for Meteorology   CDO       @Climate Data Operators version 1.6.9 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X           
�   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y           
�   time               standard_name         time   bounds        	time_bnds      units         hours since 1850-01-31 23:52:00    calendar      proleptic_gregorian         
�   	time_bnds                            
�   tsurf                         	long_name         surface temperature    units         K      code         �   table            �   
_FillValue        ����   missing_value         ����        
�                @�             @�P     C���@��     @�     @�b     C��q@�     @�     @��     C�ť@ݑ     @ٰ     @��     C��@��    @�     @�	     C��s@�V�    @�f     @�P     C��@렀    @�     @�     C���@��    @��     @���    C���@�@    @�     @�     C���@�:�    @�B�    @�7�    C��
@�_�    @�f     @�\�    C��S@��@    @��     @��     C��9@���    @���    @���    C��,@��@    @��     @��     C��	@��@    @���    @��     C���A �`    A @    A�    C���A�     A     A�    C��A��    A0�    A+@    C��A�`    AB�    A=�    C��0A�     AU     AO�    C���A��    Af�    Aa@    C��	A��    Ax�    As     C��wA     A�@    A��    C���A	�    A��    A	�@    C���A
*�    A	��    A
�     C���A<`    A
�@    A��    C��AN�    A�     A�@    C��CA`�    A�    A�     C���Ar`    A�@    A��    C���A�     A     A@    C��QAKP    A�    A��    C�˒A�0    A�     A`    C���A]    A     A�@    C��A��    A��    A%     C��Ao0    A0�    A�`    C���A�    A�     A7@    C��A��    AB�    A�     C���A	�    A��    AI     C���A�    AT�    A�@    C���A�    A��    A[     C��ZA��    Af�    A�     C��?A-�    A�    Al�    C���A��    Ax�    A�     C��)A?�    A�    A     C�ʕAȰ    A��    A�    C�ܳAQ�    A�    A��    C�թA��    A�`    A     C��Ac�    A%�    A��    C��,A�    A��    A+�    C���Aup    A7`    A��    C���A�P    A�@    A=�    C�αA�0    AI     A�`    C���A    A�     AO@    C���A��    AZ�    A�     C���A"0    A��    Aa`    C��}A�    Am     A�@    C��GA3�    A��    As     C�լA��    A~�    A�     C��HAF    A�    A�@    C��A��    A��    A     C���A +�    A �    A K�    C��A pX    A QP    A ��    C�fA ��    A ��    A Ԑ    C��A �h    A �`    A!     C�A!=�    A!�    A!]p    C��lA!�H    A!c@    A!��    C��EA!��    A!��    A!�    C��A"X    A!�P    A"*�    C��A"O�    A"0�    A"o`    C��HA"�8    A"u0    A"��    C�� A"��    A"��    A"�p    C��~A#H    A"�@    A#<�    C���A#a�    A#B�    A#�P    C��eA#�(    A#�     A#��    C�ߗA#��    A#ː    A$
`    C��SA$/8    A$0    A$N�    C��ZA$s�    A$T�    A$�@    C�؅A$�    A$�    A$װ    C��\A$��    A$݀    A%P    C��mA%A(    A%"     A%`�    C��|A%��    A%f�    A%�0    C��xA%�    A%�     A%�    C���A&�    A%�p    A&.@    C�ŅA&S    A&4    A&r�    C���A&��    A&x�    A&�     C���A&��    A&��    A&��    C��wA' �    A'`    A'@0    C��A'e    A'F     A'��    C��A'�x    A'�p    A'�    C�׶A'��    A'��    A(�    C��A(2�    A(P    A(R     C��A(v�    A(W�    A(��    C��]A(�h    A(�`    A(�     C��A(��    A(��    A)p    C��*A)Dx    A)%@    A)d    C��IA)��    A)i�    A)��    C��A)�X    A)�P    A)��    C��A*�    A)��    A*1`    C��A*Vh    A*70    A*v     C��`A*��    A*{�    A*�p    C��bA*�H    A*�@    A*��    C��hA+#�    A+�    A+CP    C��EA+hX    A+I     A+��    C��~A+��    A+��    A+�`    C��A+�8    A+�0    A,�    C���A,5�    A,�    A,U@    C���A,zH    A,[    A,��    C���A,��    A,��    A,�P    C��SA-(    A,�     A-"�    C���A-G�    A-(�    A-g0    C��A-�8    A-m     A-��    C� OA-Ш    A-��    A-�@    C�RA.    A-�    A.4�    C���A.Y�    A.:�    A.y     C��&A.�(    A.~�    A.��    C�� A.�    A.Ð    A/0    C��A/'    A/     A/F�    C���A/kx    A/Lp    A/�    C��fA/�    A/��    A/ϰ    C�A/�    A/Հ    A0
    C��!A0|    A0�    A0,H    C��A0>�    A0/0    A0N�    C���A0a    A0Qh    A0p�    C��}A0�<    A0s�    A0�    C���A0�t    A0��    A0�@    C��A0Ǭ    A0�(    A0�x    C�gA0��    A0�`    A0��    C���A14    A0��    A1     C�EA1.l    A1�    A1>8    C��A1P�    A1A     A1`p    C�aA1r�    A1cX    A1��    C��A1�,    A1��    A1��    C��A1�d    A1��    A1�0    C�KA1ٜ    A1�    A1�h    C��A1��    A1�P    A2�    C��A2$    A2�    A2-�    C��A2@\    A20�    A2P(    C�hA2b�    A2S    A2r`    C� {A2��    A2uH    A2��    C�~A2�    A2��    A2��    C�%�A2�T    A2��    A2�     C�-jA2�    A2�    A2�X    C��A3�    A2�@    A3�    C��-A30    A3 �    A3?�    C��A3RL    A3B�    A3b    C��A3t�    A3e     A3�P    C� ;A3��    A3�8    A3��    C�$tA3�    A3��    A3��    C�"0A3�D    A3��    A3�    C�XA3�|    A3��    A4H    C�'VA4�    A40    A4/�    C�"LA4B    A42�    A4Q�    C�%:A4d<    A4T�    A4t    C�~A4�t    A4v�    A4�@    C�0A4��    A4�(    A4��    C�7,A4��    A4�x    A4��    C�;�A4�4    A4ݰ    A4�     C�LJA5l    A4��    A58    C�M�A51�    A5"     A5A�    C�&�A5S�    A5Dp    A5c�    C��A5v,    A5f�    A5��    C�+A5�d    A5��    A5�0    C�.�A5��    A5�    A5ʀ    C�<�A5��    A5�h    A5�    C�.�A5�$    A5�    A6�    C�O�A6!\    A6�    A61(    C�\�