"""Build a synthetic knee CT (DICOM, LPS, left knee) from the MRI bone meshes of Dryad specimen S193761.
Bone = 1300 HU cortex / 300 HU cancellous, soft tissue 40 HU, air -1000 HU, Gaussian noise 25 HU,
plus one small detached fragment near the lateral tibial rim to exercise fragment detection.
Requires: numpy, scipy, trimesh, pydicom."""
import numpy as np,trimesh,pickle,json,os,zipfile,datetime
import pydicom
from pydicom.dataset import Dataset,FileDataset,FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian,generate_uid,CTImageStorage
from scipy import ndimage
# Usage: python make_synthetic_ct.py /path/to/STLs   (the STLs.zip folder from Dryad doi:10.5061/dryad.zcrjdfnpv)
import sys
STL_DIR=sys.argv[1] if len(sys.argv)>1 else "STLs"
D=STL_DIR+"/Raw MRI Scan STLs/S193761_Left_";O=np.array([134.17,64.73,-476.8])   # O: tibial plateau centre used as origin
def L(k):
    m=trimesh.load(D+k+"_smooth.stl"); v=m.vertices-O; v=np.c_[v[:,0],-v[:,1],v[:,2]]; return trimesh.Trimesh(v,m.faces[:,::-1])
bones={k:L(k) for k in ['Bone_Femur','Bone_Tibia','Bone_Fibula','Bone_Patella']}
# template frame -> LPS (left knee): x_L = x_t, y_P = -y_t, z_S = z_t ; then shift so values positive-ish
def to_lps(v): return np.c_[v[:,0],-v[:,1],v[:,2]]
px,sl=0.7,1.0
xs=np.arange(-80,90,px); ys=np.arange(-80,80,px); zs=np.arange(-140,150,sl)
vol=np.full((len(zs),len(ys),len(xs)),-1000,np.float32)
X,Y=np.meshgrid(xs,ys)
body=((X-5)/78)**2+((Y)/70)**2<1
vol[:,body]=40
mask_all=np.zeros(vol.shape,bool)
for k,m in bones.items():
    mv=trimesh.Trimesh(to_lps(m.vertices),m.faces[:,::-1])
    vox=mv.voxelized(pitch=0.5).fill()
    pts=vox.points  # filled voxel centers
    ix=np.round((pts[:,0]-xs[0])/px).astype(int); iy=np.round((pts[:,1]-ys[0])/px).astype(int); iz=np.round((pts[:,2]-zs[0])/sl).astype(int)
    ok=(ix>=0)&(ix<len(xs))&(iy>=0)&(iy<len(ys))&(iz>=0)&(iz<len(zs))
    mk=np.zeros(vol.shape,bool); mk[iz[ok],iy[ok],ix[ok]]=True; mk=ndimage.binary_closing(mk,iterations=1)
    mask_all|=mk; print(k,mk.sum())
# cancellous 300, cortical shell 1300
er=ndimage.binary_erosion(mask_all,iterations=2)
vol[mask_all]=1300; vol[er]=300
# synthetic Segond fragment: small ellipsoid just off the lateral tibial rim (template ALL insertion ~ (39,4,-8))
tb=bones['Bone_Tibia']; tl=trimesh.Trimesh(to_lps(tb.vertices),tb.faces[:,::-1])
guess=np.array([40.0,2.0,-8.0]); cp,_,tid=trimesh.proximity.closest_point(tl,[guess]); nrm=tl.face_normals[tid[0]]
c=cp[0]+nrm*5.0; print('fragment centre',c.round(1),'normal',nrm.round(2))
Z,Yg,Xg=np.meshgrid(zs,ys,xs,indexing='ij')
frag=((Xg-c[0])/2.2)**2+((Yg-c[1])/4.5)**2+((Zg:=Z-c[2])/3.5)**2<1
vol[frag]=1100
vol+=np.random.default_rng(0).normal(0,25,vol.shape).astype(np.float32)
print('vol',vol.shape)
out='synthetic_knee_CT'; os.makedirs(out,exist_ok=True)
for f in os.listdir(out): os.remove(os.path.join(out,f))
study,series,frame=generate_uid(),generate_uid(),generate_uid()
for i,z in enumerate(zs):
    meta=FileMetaDataset(); meta.MediaStorageSOPClassUID=CTImageStorage; meta.MediaStorageSOPInstanceUID=generate_uid(); meta.TransferSyntaxUID=ExplicitVRLittleEndian
    ds=FileDataset(None,{},file_meta=meta,preamble=b"\0"*128)
    ds.SOPClassUID=CTImageStorage; ds.SOPInstanceUID=meta.MediaStorageSOPInstanceUID; ds.Modality='CT'; ds.PatientName='SYNTHETIC^KNEE'; ds.PatientID='SYNTH001'
    ds.StudyInstanceUID=study; ds.SeriesInstanceUID=series; ds.FrameOfReferenceUID=frame; ds.SeriesDescription='Synthetic knee CT 0.7 mm bone (with Segond-like fragment)'
    ds.ImagePositionPatient=[float(xs[0]),float(ys[0]),float(z)]; ds.ImageOrientationPatient=[1,0,0,0,1,0]; ds.PixelSpacing=[px,px]; ds.SliceThickness=sl
    ds.Rows,ds.Columns=vol.shape[1],vol.shape[2]; ds.SamplesPerPixel=1; ds.PhotometricInterpretation='MONOCHROME2'; ds.BitsAllocated=16; ds.BitsStored=16; ds.HighBit=15; ds.PixelRepresentation=1
    ds.RescaleIntercept=0; ds.RescaleSlope=1; ds.Laterality='L'; ds.InstanceNumber=i+1; ds.BodyPartExamined='KNEE'; ds.ConvolutionKernel='BONE'
    ds.PixelData=np.clip(vol[i],-1024,3000).astype(np.int16).tobytes()
    ds.is_little_endian=True; ds.is_implicit_VR=False
    ds.save_as(os.path.join(out,f'IM{i:04d}.dcm'),write_like_original=False)
zf='synthetic_knee_CT.zip'
with zipfile.ZipFile(zf,'w',zipfile.ZIP_DEFLATED) as z:
    for f in sorted(os.listdir(out)): z.write(os.path.join(out,f),'synthetic_knee_CT/'+f)
print(os.path.getsize(zf)/1e6,'MB')
json.dump({"frag_center_lps":c.tolist()},open('synthetic_knee_CT_meta.json','w'))
