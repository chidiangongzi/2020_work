load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;begin

 a  = addfile("./wrfout_d02_2011-05-05_00_00_00","r")
 wks = gsn_open_wks("png","./wmax/W-MAX")
 cmap=(/(/255,255,255/),(/0,0,0/),\
       (/255,255,253/),(/255,255,253/),\
       (/255,240,187/),(/255,172,117/),\
       (/255,120,86/),(/255,61,61/),\
       (/246,39,53/),(/216,21,47/),\
       (/166,0,32/)/)/255.0
gsn_define_colormap(wks,cmap) 

times=wrf_user_getvar(a,"times",-1)
ntimes=dimsizes(times)

;qcloud=a->QCLOUD(:,0:28,99:209,55:200)
;qrain=a->QRAIN(:,0:28,99:209,55:200)
lat0              := a->XLAT
lat:=lat0(0,:,0)

lon0            = a->XLONG
lon=lon0(0,0,:)

w0=wrf_user_getvar(a, "wa", -1)

w_00=new((/49,29,468,471/), "float")
w_01=new((/49,29,468,471/), "float")

 p1 = wrf_user_getvar(a,"pressure",-1)
  ua1 = wrf_user_getvar(a,"ua",-1)
  va1 = wrf_user_getvar(a,"va",-1)
var2d_u1=wrf_interp_3d_z(ua1, p1, 650.)
var2d_v1=wrf_interp_3d_z(va1, p1, 650.)
;var2d_u1=ua1(:,3,:,:)
;var2d_v1=va1(:,3,:,:)
var2d_u1@_FillValue=-999
var2d_v1@_FillValue=-999

do i = 0, 470,1
    do j = 0, 467,1
        do k = 0, 28,1
            do m = 0, 48,1

             if (w0(m,k,j,i) .gt. 0) then
                w_00(m,k,j,i)=w0(m,k,j,i)
            else
                w_00(m,k,j,i)=-999
            
             end if
                ; write loop content
            end do
            ; write loop content
        end do
        ; write loop content
    end do
    ; write loop content
end do

w_00@_FillValue=-999
w_00!2="south_north"
w_00&south_north=lat
w_00!3="west_east"
w_00&west_east=lon

 nw1=dimsizes(w_00)
  ntimes=nw1(0)
  nlatis=nw1(2)
  nlongis=nw1(3)
  max_w2=new((/ntimes,nlatis,nlongis/), "float")

do t= 0, ntimes-1,1
    do i2= 0, nlatis-1,1
      do j2 = 0, nlongis-1,1
        max_w2(t,i2,j2)=max(w_00(t,:,i2,j2))
      end do
    end do
  end do

max_w2@_FillValue=-999
max_w2!1="south_north"
max_w2&south_north=lat
max_w2!2="west_east"
max_w2&west_east=lon

lat2d             = a->XLAT    
lon2d             = a->XLONG

cnres                   =  True              ; contour/map resources
  vcres                   =  True              ; vector resources

  cnres@gsnDraw           = False              ; Turn these off. We
  cnres@gsnFrame          = False              ; will overlay plots
  vcres@gsnDraw           = False              ; later.
  vcres@gsnFrame          = False
;
; Lambert conformal projections are limited using
; the corners method rather than the latlon method
; seen for cylindrical equidistant projections.
;
   cnres@mpOutlineOn                 = True    
cnres@mpDataBaseVersion       = "MediumRes" 
cnres@mpDataSetName           = "Earth..4"  
cnres@mpProvincialLineThicknessF =4 

cnres@mpGeophysicalLineThicknessF= 2.              
cnres@mpNationalLineThicknessF= 2.                 
cnres@mpOutlineBoundarySets = "AllBoundaries"

cnres@pmTickMarkDisplayMode = "Always"


    cnres@mpProjection      = "CylindricalEquidistant"   
     cnres@mpMinLatF         = 40
     cnres@mpMaxLatF         = 47.3
     cnres@mpMinLonF         = 120.5
     cnres@mpMaxLonF         = 132
     ;cnres@mpMinLatF         = 27.5
     ;cnres@mpMaxLatF         = 32
     ;cnres@mpMinLonF         = 117
     ;cnres@mpMaxLonF         = 122
     cnres@mpCenterLonF      = (cnres@mpMinLonF + cnres@mpMaxLonF) / 2.
cnres@cnLevelSelectionMode="ManualLevels"

  cnres@cnMinLevelValF=0
  cnres@cnMaxLevelValF=8
  cnres@cnLevelSpacingF=1
  cnres@gsnAddCyclic            = False            ; regional data 
  cnres@cnFillOn                = True
  cnres@cnLinesOn               = False           ; turn off contour lines


  vcres@vcRefMagnitudeF          = 10.0             ; define vector ref mag
  vcres@vcRefLengthF             = 0.045            ; define length of vec ref
  vcres@vcGlyphStyle             = "CurlyVector"    ; turn on curly vectors
  vcres@vcMinDistanceF           = 0.017            ; thin vectors
  vcres@vcRefAnnoString1On=False
  vcres@vcRefAnnoString2On=True
  vcres@vcRefAnnoString2="10 m/s"
  vcres@vcRefAnnoSide="Top"
  vcres@vcRefAnnoOrthogonalPosF=-0.12
  vcres@vcRefAnnoParallelPosF=0.95
  ;vcres@vcRefAnnoOrthogonalPosF  = -1.6               ; move ref vector down
  vcres@gsnAddCyclic             = False            ; regional data 

;---Make sure vector plot doesn't have subtitles
  vcres@gsnLeftString     = ""
  vcres@gsnRightString    = ""


do nt = 0,ntimes-1,1 
     cnres@tiMainString         = times(nt)


     var2d_u11=var2d_u1(nt,:,:)
     var2d_v11=var2d_v1(nt,:,:)
     var2d_u11@lat2d             = a->XLAT(nt,:,:)    
     var2d_u11@lon2d             = a->XLONG(nt,:,:)
     var2d_v11@lat2d             = a->XLAT(nt,:,:)    
     var2d_v11@lon2d             = a->XLONG(nt,:,:)
     rainnc=max_w2(nt,:,:)
     rainnc@lat2d             = a->XLAT(nt,:,:)    
     rainnc@lon2d             = a->XLONG(nt,:,:)

     contour_fill_plot = gsn_csm_contour_map(wks,rainnc,cnres)
     vector_plot       = gsn_csm_vector(wks,var2d_u11,var2d_v11,vcres)
     ;plot = gsn_csm_contour(wks,rainnc,resc)

   overlay(contour_fill_plot,vector_plot)

  draw(contour_fill_plot)    ; This will draw everything
  frame(wks)


end do
