program Trajectory_into_pdb
    implicit none

    integer :: natm,i,n,k,l,p,NSS,m                  !!!! NSS number of snapshots
    integer, dimension(:),allocatable :: ASN         !!!!ASN= Atom Serial Number (Second column of PDB)
    character(4), dimension(:),allocatable :: ATOM   !! First column of the PDB
    character(4), dimension(:),allocatable :: ATOM_Name   !! Third column of the PDB
    character(3), dimension(:),allocatable :: Res_Name   !! 4th column of the PDB
    character(1), dimension(:),allocatable :: Ci   !! Ci=Chain Identifier (5th column of the PDB)
    integer, dimension(:),allocatable :: RSN   !! RSN=Residue sequence number (6th column of the PDB)
    real, dimension(:),allocatable :: xcor1, ycor1,zcor1 !coordinates xyz (7th, 8th and 9th column of the PDB)
    real, dimension(:),allocatable :: N1,N2  !(10th and 11th column of the PDB)
    character(1), dimension(:),allocatable :: com   !com=component (10th and 11th column of the PDB)
    character(10), dimension(:),allocatable :: atomT,row1T, AidT
    real, dimension(:),allocatable ::  XT,YT,ZT
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    ! using keyword arguments
    call date_and_time(date,time,zone,values)
    call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)
    call date_and_time(VALUES=values)
    write(*,105) date, time, zone

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call length(natm)
    allocate(ASN(natm),ATOM(natm),ATOM_Name(natm),Res_Name(natm),Ci(natm),RSN(natm),xcor1(natm),ycor1(natm),zcor1(natm))!,
    allocate(N1(natm),N2(natm),com(natm),atomT(natm),row1T(natm),AidT(natm),XT(natm),YT(natm),ZT(natm))

    open(unit=1, file= "Trajectory_input.txt")
    open(unit=3, file= "pdb_SAMPLE.pdb")
    open(unit=4, file= "PDB_file.pdb")
    write(*,*) natm

    NSS=2001  !Number of snapshots must be manually changed


 do p=1,natm

 read(3,*) ATOM(p),ASN(p),ATOM_Name(p), Res_Name(p), Ci(p),RSN(p),xcor1(p),ycor1(p),zcor1(p),N1(p),N2(p),com(p)

end do
close(3)
write(4,102) "HEADER    2021-02-17 13:33:00.346718"

read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
do i=1,NSS
 write(*,106) i
  write(4,106) i
     Do k=1,natm

read(1,*) row1T(k),AidT(k),atomT(k), XT(k),YT(k),ZT(k)

!!!!!!!!!!!!!! If- condition definitions!!!!!!!!!!!!!!!
    if (ATOM_Name(k) .EQ. "N") THEN
    ATOM_Name(k)=" N"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ATOM_Name(k) .EQ. "CA") THEN
    ATOM_Name(k)=" CA"
    END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ATOM_Name(k) .EQ. "OT1") THEN
    ATOM_Name(k)=" OT1"
    END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ATOM_Name(k) .EQ. "OT2") THEN
    ATOM_Name(k)=" OT2"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (ATOM_Name(k) .EQ. "C") THEN
    ATOM_Name(k)=" C"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (ATOM_Name(k) .EQ. "O") THEN
    ATOM_Name(k)=" O"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "H1") THEN
    ATOM_Name(k)=" H1"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "H2") THEN
    ATOM_Name(k)=" H2"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "H3") THEN
    ATOM_Name(k)=" H3"
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HA2") THEN
    ATOM_Name(k)=" HA2"
    END IF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HA1") THEN
    ATOM_Name(k)=" HA1"
    END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HA3") THEN
    ATOM_Name(k)=" HA3"
    END IF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CB") THEN
    ATOM_Name(k)=" CB"
    END IF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CG") THEN
    ATOM_Name(k)=" CG"
    END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OH") THEN
    ATOM_Name(k)=" OH"
    END IF
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HH") THEN
    ATOM_Name(k)=" HH"
    END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CG1") THEN
    ATOM_Name(k)=" CG1"
    END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CD") THEN
    ATOM_Name(k)=" CD"
    END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CE") THEN
    ATOM_Name(k)=" CE"
    END IF
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "NZ") THEN
    ATOM_Name(k)=" NZ"
    END IF
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "H") THEN
    ATOM_Name(k)=" H"
    END IF
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HA") THEN
    ATOM_Name(k)=" HA"
    END IF
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HB1") THEN
    ATOM_Name(k)=" HB1"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HT1") THEN
    ATOM_Name(k)=" HT1"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HT2") THEN
    ATOM_Name(k)=" HT2"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HT3") THEN
    ATOM_Name(k)=" HT3"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HB2") THEN
    ATOM_Name(k)=" HB2"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HB3") THEN
    ATOM_Name(k)=" HB3"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HG2") THEN
    ATOM_Name(k)=" HG2"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HG3") THEN
    ATOM_Name(k)=" HG3"
    END IF
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HG") THEN
    ATOM_Name(k)=" HG"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HD2") THEN
    ATOM_Name(k)=" HD2"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HD3") THEN
    ATOM_Name(k)=" HD3"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HE2") THEN
    ATOM_Name(k)=" HE2"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HE3") THEN
    ATOM_Name(k)=" HE3"
    END IF
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HZ") THEN
    ATOM_Name(k)=" HZ"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HZ1") THEN
    ATOM_Name(k)=" HZ1"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HZ2") THEN
    ATOM_Name(k)=" HZ2"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HZ3") THEN
    ATOM_Name(k)=" HZ3"
    END IF
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HB") THEN
    ATOM_Name(k)=" HB"
    END IF
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 if (ATOM_Name(k) .EQ. "HN") THEN
    ATOM_Name(k)=" HN"
    END IF
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CD1") THEN
    ATOM_Name(k)=" CD1"
    END IF
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CG2") THEN
    ATOM_Name(k)=" CG2"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OE1") THEN
    ATOM_Name(k)=" OE1"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    if (ATOM_Name(k) .EQ. "OE2") THEN
    ATOM_Name(k)=" OE2"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "NE2") THEN
    ATOM_Name(k)=" NE2"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CD2") THEN
    ATOM_Name(k)=" CD2"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "NE") THEN
    ATOM_Name(k)=" NE"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "CZ") THEN
    ATOM_Name(k)=" CZ"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "NH1") THEN
    ATOM_Name(k)=" NH1"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "NH2") THEN
    ATOM_Name(k)=" NH2"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "HE") THEN
    ATOM_Name(k)=" HE"
    END IF
                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OD1") THEN
    ATOM_Name(k)=" OD1"
    END IF
                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OD2") THEN
    ATOM_Name(k)=" OD2"
    END IF

                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "ND2") THEN
    ATOM_Name(k)=" ND2"
    END IF
                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OG1") THEN
    ATOM_Name(k)=" OG1"
    END IF
                                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "SD") THEN
    ATOM_Name(k)=" SD"
    END IF
                                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (ATOM_Name(k) .EQ. "OG") THEN
    ATOM_Name(k)=" OG"
    END IF
                                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "HG1") THEN
    ATOM_Name(k)=" HG1"
    END IF
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "ND1") THEN
    ATOM_Name(k)=" ND1"
    END IF
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "CE1") THEN
    ATOM_Name(k)=" CE1"
    END IF
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "CE2") THEN
    ATOM_Name(k)=" CE2"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "CZ") THEN
    ATOM_Name(k)=" CZ"
    END IF
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "HD1") THEN
    ATOM_Name(k)=" HD1"
    END IF
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "HE1") THEN
    ATOM_Name(k)=" HE1"
    END IF
                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (ATOM_Name(k) .EQ. "OXT") THEN
    ATOM_Name(k)=" OXT"
    END IF
 write(4,101) ASN(k),ATOM_Name(k), Res_Name(k), Ci(k),RSN(k),XT(k), YT(k),ZT(k),N1(k),N2(k),ATOM_Name(k)
 write(*,101) ASN(k),ATOM_Name(k), Res_Name(k), Ci(k),RSN(k),XT(k), YT(k),ZT(k),N1(k),N2(k),ATOM_Name(k)

end do
         if (i .NE. NSS) THEN


read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
    END IF
          write(*,102) "TER     139      GLY A  10"
          write(4,102) "TER     139      GLY A  10"
          write(4,102) "ENDMDL"
          write(*,102) "ENDMDL"
            if (i .eq. NSS) THEN
            write(*,102)"MASTER      102    0    0    0    2    0    0    6   77    1    0    1"
            write(4,102)"MASTER      102    0    0    0    2    0    0    6   77    1    0    1"
            write(*,102)"END"
            write(4,102)"END"
            end if



end do

    pause
100 format(A4,2x,A4,2x,A4,2x, F7.3,2x, F8.3,2x, F8.3)
101 format("ATOM  ",I5,1X,A4,x,A3,X,A1,I4,4X,3F8.3,2F6.2,10X,A2,A2)
102 format(a)
105 format("Generated by Amirhossein Saali    ",a,2x,a,2x,a,'                 ')
106 format("MODEL        ",I6)
end program Trajectory_into_pdb
