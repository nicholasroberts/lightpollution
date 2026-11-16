"""Extract linear camera-blue Bayer samples; no rendered RGB or white balance.

Default exposure mapping is museum DSC_5245..DSC_5256, verified against the
uploaded TIFF XMP. Other datasets require an explicit mapping JSON. This is
a preparation step, not registration or a complete instrument calibration.
Dependencies: numpy, rawpy, exifread, tifffile, opencv-python-headless.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

MUSEUM_TIMES=['22:09:51','22:10:04','22:10:17','22:11:30','22:11:43','22:11:56',
              '22:12:29','22:12:41','22:12:54','22:13:44','22:13:57','22:14:10']
MUSEUM_MAPPING=[dict(exposure=n,stem=f'DSC_{5244+n}',
    datetime_original='2019:05:14 '+MUSEUM_TIMES[n-1]) for n in range(1,13)]


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def select_sources(folder,mapping,raw_kind='auto'):
    files=[p for p in Path(folder).iterdir() if p.is_file() and p.suffix.lower() in ('.nef','.dng')]
    result=[]
    for record in mapping:
        stem=record['stem'].casefold()
        kinds=['.nef','.dng'] if raw_kind=='auto' else ['.'+raw_kind]
        found=None
        for suffix in kinds:
            matches=[p for p in files if p.stem.casefold()==stem and p.suffix.lower()==suffix]
            if len(matches)>1:raise ValueError('Ambiguous source: '+record['stem']+suffix)
            if matches:found=matches[0];break
        if found is None:raise ValueError('Missing original RAW: '+record['stem']+' (NEF or DNG)')
        result.append((record,found))
    if len({p.resolve() for _,p in result})!=len(result):raise ValueError('RAW source used more than once')
    return result


def _number(tag):
    if tag is None:raise ValueError('Required RAW exposure metadata absent')
    value=tag.values[0]
    if hasattr(value,'num') and hasattr(value,'den'):return float(value.num)/float(value.den)
    return float(value)


def read_metadata(path):
    import exifread
    with Path(path).open('rb') as stream:tags=exifread.process_file(stream,details=False)
    time=tags.get('EXIF DateTimeOriginal')
    if time is None:raise ValueError('RAW DateTimeOriginal absent: '+str(path))
    return dict(datetime_original=str(time),
        exposure_seconds=_number(tags.get('EXIF ExposureTime')),
        iso=_number(tags.get('EXIF ISOSpeedRatings') or tags.get('EXIF PhotographicSensitivity')),
        f_number=_number(tags.get('EXIF FNumber')),
        model=str(tags.get('Image Model','unknown')),
        subsecond_original=str(tags.get('EXIF SubSecTimeOriginal','')))


def blue_bayer_plane(mosaic,colors,black_levels,color_desc):
    """One real blue sample per Bayer 2x2 cell; no interpolation or gain.

    Colors must correspond to the *visible* mosaic, including its CFA phase.
    Negative black-subtracted values are retained, not clipped or normalized.
    """
    mosaic=np.asarray(mosaic);colors=np.asarray(colors)
    if mosaic.ndim!=2 or colors.shape!=(2,2):raise ValueError('Expected a visible 2D Bayer mosaic and its 2x2 CFA phase')
    desc=color_desc.decode('ascii') if isinstance(color_desc,bytes) else str(color_desc)
    positions=[(y,x,int(colors[y,x])) for y in range(2) for x in range(2)
               if int(colors[y,x])<len(desc) and desc[int(colors[y,x])]=='B']
    if len(positions)!=1:raise ValueError('Expected exactly one blue site per 2x2 Bayer cell')
    oy,ox,index=positions[0];h,w=(mosaic.shape[0]//2,mosaic.shape[1]//2)
    if h<2 or w<2:raise ValueError('RAW mosaic too small')
    plane=mosaic[oy:2*h:2,ox:2*w:2].astype(np.float32)-float(black_levels[index])
    return plane,dict(blue_color_index=index,blue_offset_yx=[oy,ox],blue_black_level=float(black_levels[index]),
                      dropped_bottom_rows=int(mosaic.shape[0]%2),dropped_right_columns=int(mosaic.shape[1]%2))


def output_shape(shape,long_edge):
    if long_edge<=0:return tuple(shape)
    ratio=min(1.,long_edge/max(shape))
    return tuple(max(1,int(round(n*ratio))) for n in shape)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-dir',type=Path,required=True,help='One folder containing the original NEF or DNG files')
    p.add_argument('--output',type=Path,help='Fresh output directory; originals are read only')
    p.add_argument('--inspect',action='store_true',help='Print RAW settings without converting or creating files')
    p.add_argument('--raw-kind',choices=['auto','nef','dng'],default='auto',help='auto prefers NEF if both formats exist')
    p.add_argument('--mapping',type=Path,help='JSON list of {exposure,stem,datetime_original}; default is the verified museum mapping')
    p.add_argument('--long-edge',type=int,default=2000,help='Common area-downsampled long edge; 0 retains the blue Bayer grid')
    a=p.parse_args()
    if not a.inspect and a.output is None:p.error('--output is required unless --inspect is used')
    if not a.inspect and a.output.exists():p.error('Use a fresh output directory')
    if not a.raw_dir.is_dir():p.error('RAW folder not found')
    if a.long_edge<0:p.error('--long-edge must be nonnegative')
    mapping=json.loads(a.mapping.read_text()) if a.mapping else MUSEUM_MAPPING
    if not isinstance(mapping,list) or sorted(r['exposure'] for r in mapping)!=list(range(1,13)):
        p.error('Mapping must contain exactly exposures 1..12')
    mapping=sorted(mapping,key=lambda r:r['exposure'])
    try:
        sources=select_sources(a.raw_dir,mapping,a.raw_kind)
        prepared=[]
        for record,path in sources:
            metadata=read_metadata(path)
            if metadata['datetime_original']!=record['datetime_original']:
                raise ValueError('RAW capture timestamp differs from mapped TIFF: '+path.name)
            prepared.append((record,path,metadata))
        print('Exposure  File                  Time                 Seconds       ISO       f-number',flush=True)
        for record,path,metadata in prepared:
            print(f"{record['exposure']:8d}  {path.name:20s}  {metadata['datetime_original']:19s}  "
                  f"{metadata['exposure_seconds']:10g}  {metadata['iso']:8g}  {metadata['f_number']:8g}",flush=True)
        if a.inspect:return
        import rawpy,cv2,tifffile
        first=prepared[0][2]
        differences=[]
        for _,path,metadata in prepared[1:]:
            changed=[f"{k}: {first[k]:g} -> {metadata[k]:g}" for k in ('iso','exposure_seconds','f_number')
                     if not np.isclose(metadata[k],first[k],rtol=0,atol=1e-8)]
            if changed:differences.append(path.name+' ('+', '.join(changed)+')')
            if metadata['model']!=first['model']:raise ValueError('Camera model changes within the set')
        if differences:raise ValueError('Exposure settings differ from '+prepared[0][1].name+': '+'; '.join(differences))
    except ImportError as exc:
        p.error('Install numpy rawpy exifread tifffile opencv-python-headless: '+str(exc))
    except (ValueError,KeyError,OSError) as exc:p.error(str(exc))
    a.output.mkdir(parents=True)
    manifest=dict(schema_version=1,status='linear_raw_blue_prepared_not_registered',
        mapping_origin='explicit_json' if a.mapping else 'museum_TIFF_XMP_DSC_5245_to_DSC_5256',
        processing=dict(method='visible_Bayer_blue_samples_minus_camera_reported_black',
            units='RAW_digital_counts_above_camera_black',white_balance_applied=False,gamma_applied=False,
            color_matrix_applied=False,display_profile_or_look_applied=False,auto_brightness_applied=False,
            demosaicing_applied=False,negative_values_retained=True,per_frame_brightness_normalization=False,
            resize='one_common_area_resampling' if a.long_edge else 'none',requested_long_edge=a.long_edge),
        analyzer_groups=[[1,4,7,10],[2,5,8,11],[3,6,9,12]],
        software=dict(rawpy=rawpy.__version__,numpy=np.__version__,opencv=cv2.__version__,tifffile=tifffile.__version__),
        script_sha256=sha256(__file__),input_shape=None,output_shape=None,frames=[],
        limitations=['Not a full polarimetric instrument calibration.',
            'Blue Bayer extraction is a single camera channel, not rendered Adobe RGB blue.',
            'New raw decoding/CFA sampling can change pixel geometry; re-register this set.',
            'Do not substitute these files into old checksum-locked registration folders.',
            'Astronomical anchors, theoretical patterns and ALAN settings are not fitted or modified.'])
    try:
        for record,path,metadata in prepared:
            with rawpy.imread(str(path)) as raw:
                mosaic=raw.raw_image_visible
                if mosaic is None:raise ValueError('No visible Bayer mosaic: '+path.name)
                colors=raw.raw_colors_visible[:2,:2].copy()
                plane,info=blue_bayer_plane(mosaic,colors,raw.black_level_per_channel,raw.color_desc)
                input_shape=list(mosaic.shape);shape=output_shape(plane.shape,a.long_edge)
                if manifest['input_shape'] is None:
                    manifest['input_shape']=input_shape;manifest['output_shape']=list(shape)
                elif input_shape!=manifest['input_shape'] or list(shape)!=manifest['output_shape']:
                    raise ValueError('RAW geometry differs within the set: '+path.name)
                if manifest['frames'] and (colors.tolist()!=manifest['frames'][0]['cfa_visible_2x2']
                    or info['blue_offset_yx']!=manifest['frames'][0]['blue_offset_yx']):
                    raise ValueError('Visible CFA phase differs within the set: '+path.name)
                oy,ox=info['blue_offset_yx'];h,w=plane.shape
                white=float(raw.white_level)
                saturation_fraction=float(np.mean(mosaic[oy:2*h:2,ox:2*w:2]>=white))
                if plane.shape!=shape:plane=cv2.resize(plane,(shape[1],shape[0]),interpolation=cv2.INTER_AREA)
                output=a.output/f"{record['exposure']}.tif"
                description=json.dumps(dict(kind='linear_raw_camera_blue',channel=None,
                    datetime_original=metadata['datetime_original'],source_file=path.name,
                    camera_blue_color_index=info['blue_color_index'],no_display_encoding=True))
                tifffile.imwrite(output,plane.astype(np.float32),photometric='minisblack',
                    compression='zlib',metadata=None,description=description)
                frame=dict(exposure=record['exposure'],source_path=str(path.resolve()),source_sha256=sha256(path),
                    output_path=str(output.resolve()),output_sha256=sha256(output),channel=None,
                    dtype='float32',raw_white_level=white,raw_blue_saturated_fraction=saturation_fraction,
                    cfa_visible_2x2=colors.tolist(),black_level_per_channel=list(raw.black_level_per_channel),
                    color_desc=raw.color_desc.decode('ascii'),**metadata,**info)
                manifest['frames'].append(frame)
                (a.output/'linear_raw_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
                print(f"{record['exposure']:2d}/12 {path.name} -> {output.name}; shape {plane.shape}",flush=True)
        manifest['status']='linear_raw_blue_prepared_not_registered_complete'
        (a.output/'linear_raw_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    except Exception as exc:
        manifest['status']='conversion_incomplete';manifest['error']=str(exc)
        (a.output/'linear_raw_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        raise
    print('Prepared 12 linear camera-blue TIFFs and linear_raw_manifest.json. Re-registration is required.')


if __name__=='__main__':main()

