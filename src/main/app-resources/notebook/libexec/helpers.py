import cioppy
import geopandas as gp
import os, shutil
import pandas as pd

from shapely.geometry import box
from snappy import jpy
from snappy import ProductIO
from snappy import GPF
import snappy
import gdal
import osr
import ogr
from shapely.geometry import box
import json
import sys


def get_metadata(input_references, data_path):
    
    ciop = cioppy.Cioppy()

    if isinstance(input_references, str):

        search_params = dict()

        search_params['do'] = 'terradue'

        products = gp.GeoDataFrame(ciop.search(end_point=input_references, 
                            params=search_params,
                            output_fields='identifier,self,wkt,startdate,enddate,enclosure,orbitDirection,track,orbitNumber', 
                            model='EOP'))

    else:    

        temp_results = []

        for index, self in enumerate(input_references):

            search_params = dict()

            search_params['do'] = 'terradue'

            temp_results.append(ciop.search(end_point=self, 
                                params=search_params,
                                output_fields='identifier,self,wkt,startdate,enddate,enclosure,orbitDirection,track,orbitNumber', 
                                model='EOP')[0])

        products = gp.GeoDataFrame(temp_results)
        
        products = products.merge(products.apply(lambda row: analyse(row, data_path), axis=1),
                                    left_index=True,
                                  right_index=True)
        
        
    return products

def analyse(row, data_path):
    
    series = dict()

    series['local_path'] = os.path.join(data_path, row['identifier'], row['identifier'] + '.SAFE')
   
    return pd.Series(series)


def group_analysis(df):
    df['ordinal_type'] = 'NaN'
    slave_date=df['startdate'].min()[:10]
    master_date=df['startdate'].max()[:10]
    for i in range(len(df)):
    
        if slave_date == df.iloc[i]['startdate'][:10]:
            df.loc[i,'ordinal_type']='Pre'
    
        elif master_date == df.iloc[i]['startdate'][:10]:
            df.loc[i,'ordinal_type']='Pst'

    return 


def cog(input_tif, output_tif):
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co TILED=YES ' \
                                                                    '-co COPY_SRC_OVERVIEWS=YES ' \
                                                                    ' -co COMPRESS=LZW'))

    ds = gdal.Open(input_tif, gdal.OF_READONLY)

    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    ds.BuildOverviews('NEAREST', [2,4,8,16,32])
    
    ds = None

    ds = gdal.Open(input_tif)
    gdal.Translate(output_tif,
                   ds, 
                   options=translate_options)
    ds = None

    os.remove('{}.ovr'.format(input_tif))
    os.remove(input_tif)
    
    
def mosaic_inputs(input_prod): 
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    
    slices = jpy.array('org.esa.snap.core.datamodel.Product', input_prod.identifier.count())
    for index, row in input_prod.iterrows(): 
        slices[index] = snappy.ProductIO.readProduct("%s/MTD_MSIL2A.xml" % row['local_path'])

    mosaic = snappy.GPF.createProduct('SliceAssembly', parameters, slices)
    return mosaic


def resample2ref_band(product,reference_band):
    
    snappy.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    parameters = HashMap()
    parameters.put('referenceBand', reference_band)
    product = snappy.GPF.createProduct('Resample', parameters, product)
    return product


def subset_to_aoi_reduce_bands(product,wkt,req_bands):
    snappy.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    parameters = HashMap()
    parameters.put('sourceBands', ','.join(req_bands))
    parameters.put('geoRegion', wkt)
    parameters.put('subSamplingX', '1')
    parameters.put('subSamplingY', '1')
    parameters.put('fullSwath', 'false')
    parameters.put('tiePointGridNames', '')
    parameters.put('copyMetadata', 'true')
    subset = snappy.GPF.createProduct('Subset', parameters, product)
    return subset

def snap_rgb(product,rgb_band_list,output_name):
    snappy.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    parameters = HashMap()
    parameters.put('sourceBands', ','.join(rgb_band_list))
    
    parameters.put('subSamplingX', '1')
    parameters.put('subSamplingY', '1')
    parameters.put('fullSwath', 'false')
    parameters.put('tiePointGridNames', '')
    parameters.put('copyMetadata', 'true')
    rgb = GPF.createProduct('Subset', parameters, product)
    ProductIO.writeProduct(rgb, 'tmp_rgb.tif', 'GeoTIFF-BigTIFF')
    cog('tmp_rgb.tif',output_name)
    return
    
## Replace the band name with anonymos band number     
def snap_mask(product,output_name):
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    BandDescriptor = snappy.jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')
    targetBand0 = BandDescriptor()
    targetBand0.name = 'mask'
    targetBand0.type = 'uint8'
    targetBand0.expression = 'quality_scene_classification<2?0:quality_scene_classification==3?0:quality_scene_classification==8?0:1'
    targetBands = snappy.jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand0
    parameters = HashMap()
    parameters.put('targetBands', targetBands)

    cloud_mask = snappy.GPF.createProduct('BandMaths', parameters,product)
    ProductIO.writeProduct(cloud_mask, 'tmp_cloud_mask.tif', 'GeoTIFF-BigTIFF')
    cog('tmp_cloud_mask.tif',output_name)
    return 
    
    
def write_tif(layer, output_name, width, height, input_geotransform, input_georef):

    driver = gdal.GetDriverByName('GTiff')

    output = driver.Create(output_name, 
                           width, 
                           height, 
                           1, 
                           gdal.GDT_Byte) 
        
    output.SetGeoTransform(input_geotransform)
    output.SetProjection(input_georef)
    output.GetRasterBand(1).WriteArray(layer)

    output.FlushCache()


def polygonize(input_tif, band, epsg, mask=None):
    
    epsg_code = int(epsg)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)

    source_raster = gdal.Open(input_tif)
    band_source = source_raster.GetRasterBand(band)

    source_mask = gdal.Open(mask)
    band_mask = source_mask.GetRasterBand(1)
    
    out_vector_file = "polygonized.json"
    driver = ogr.GetDriverByName('GeoJSON')

    out_data_source = driver.CreateDataSource(out_vector_file)

    out_layer = out_data_source.CreateLayer(out_vector_file, srs=srs)
    
    new_field = ogr.FieldDefn('change_detection', ogr.OFTInteger)
    out_layer.CreateField(new_field)

    gdal.Polygonize(band_source, band_mask, out_layer, 0, [], callback=None)

    out_data_source = None
    source_mask = None
    source_raster = None
    

    data = json.loads(open(out_vector_file).read())
    gdf = gp.GeoDataFrame.from_features(data['features'])

    gdf.crs = {'init':'epsg:{}'.format(epsg)}
    gdf = gdf.to_crs(epsg=epsg)
    
    #os.remove(out_vector_file)
    
    return gdf


def sieve_filter(input_file, output_file, pp_threshold):
    
    shutil.copy(input_file,output_file)
    ds = gdal.Open(input_file)
    output = gdal.Open(output_file,gdal.GA_Update)
    
    gdal.SieveFilter(ds.GetRasterBand(1), None, output.GetRasterBand(1), pp_threshold, 8)
    
    output.FlushCache()
    ds.FlushCache()
    ds = None
    ouput = None
    
