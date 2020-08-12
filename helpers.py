import shutil
import gdal
import os, sys, osr, ogr
import geopandas as gp
import json


def cog(ds,output_name):
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co TILED=YES ' \
                                                                    '-co COPY_SRC_OVERVIEWS=YES ' \
                                                                    '-co COMPRESS=LZW '))
        
    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    ds.BuildOverviews('NEAREST', [2,4,8,16,32])
    gdal.Translate(output_name,ds,options=translate_options)
    ds = None


def write_tif(layer, output_tif, width, height, input_geotransform, input_georef, to_cog=True):

    driver = gdal.GetDriverByName('GTiff')

    temp_out = 'temp_{}'.format(output_tif)
    
    ds = driver.Create(temp_out, 
                       width, 
                       height, 
                       1, 
                       gdal.GDT_Byte) 
        
    ds.SetGeoTransform(input_geotransform)
    ds.SetProjection(input_georef)
    ds.GetRasterBand(1).WriteArray(layer)
    
    if to_cog:
        cog(ds,output_tif)
        os.remove(temp_out)

    else:
        shutil.move(temp_out,output_tif)  
        
    ds.FlushCache()
    ds = None
    
    
    
def write_RGB(layers, output_tif, width, height, input_geotransform, input_georef, to_cog=True):

    driver = gdal.GetDriverByName('GTiff')
    
    temp_out = 'temp_{}'.format(output_tif)
    
    ds = driver.Create(temp_out, 
                       width,
                       height,
                       3, 
                       gdal.GDT_UInt16)
        
    ds.SetGeoTransform(input_geotransform)
    ds.SetProjection(input_georef)
    
    for i,layer in enumerate(layers):

        ds.GetRasterBand(i+1).WriteArray(layer)

    if to_cog:
        cog(ds,output_tif)
        os.remove(temp_out)
    
    else:
        shutil.move(temp_out,output_tif)  
        
    ds.FlushCache()
    ds = None
    
    
def polygonize(input_tif, poligonized_file, band, epsg, mask):
    
    epsg_code = epsg
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg_code))

    source_raster = gdal.Open(input_tif)
    band_raster = source_raster.GetRasterBand(band)
    
    source_mask = gdal.Open(mask)
    band_mask = source_mask.GetRasterBand(1)


    out_vector_file = poligonized_file

    driver = ogr.GetDriverByName('GeoJSON')

    out_data_source = driver.CreateDataSource(out_vector_file)
    out_layer = out_data_source.CreateLayer(out_vector_file, srs=srs)

    new_field = ogr.FieldDefn('change_detection', ogr.OFTInteger)
    out_layer.CreateField(new_field)

    gdal.Polygonize(band_raster, band_mask, out_layer, 0, [], callback=None )

    out_data_source = None
    source_raster = None
    source_mask = None
    
    
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
    