import SimpleITK as sitk

warp_f = sitk.ReadImage(snakemake.input.warp, sitk.sitkVectorFloat64)
warp_i = sitk.ReadImage(snakemake.input.invwarp, sitk.sitkVectorFloat64)
affine = sitk.ReadTransform(snakemake.input.affine)
warp = sitk.DisplacementFieldTransform(warp_f)
warp.SetInverseDisplacementField(warp_i)
composite = sitk.CompositeTransform([affine, warp])
inv_composite = sitk.CompositeTransform([warp.GetInverse(), affine.GetInverse()])
inv_composite.FlattenTransform()
composite.FlattenTransform()
inv_composite.WriteTransform(snakemake.output.inverse)
composite.WriteTransform(snakemake.output.forward)

