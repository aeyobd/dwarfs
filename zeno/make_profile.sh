
# halogsp				 Generate profile for halo model
#  out=???                         Output file for profile tables
#  m_a=1.0                         Mass within radius a
#  a=1.0                           Radial scale of model
#  b=4.0                           Radius to begin taper
#  taper=exp                       Tapering: exp, gauss, or sw
#  npoint=257                      Number of points in tables
#  rrange=1/4096:16                Range of radii tabulated
#  VERSION=1.3                     Josh Barnes  31 January 2015

if [ $# -ne 1 ]; then
  echo "Usage: $0 <filename>"
  exit 1
fi


N=1024

filename=$1

halo_path=profiles/$filename.gsp
txt_path=tmp/$filename.txt
csv_path=profiles/$filename.csv

rm -f $halo_path
#$ZENOPATH/bin/halogsp $halo_path m_a=1.0 a=1.0 npoint=$N rrange=1/16563:1024
$ZENOPATH/bin/halogsp $halo_path #rrange=1/4096:512

$ZENOPATH/bin/tsf $halo_path maxline=$N > $txt_path
julia parse_profile.jl $txt_path $csv_path
