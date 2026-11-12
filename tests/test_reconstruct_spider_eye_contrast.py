"""Scientific and file-I/O checks using synthetic skies, never field validation."""
import contextlib
import importlib.util
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import reconstruct_spider_eye_contrast as model


class EyeModelTests(unittest.TestCase):
    def setUp(self):
        self.i = np.ones((11, 11))
        self.mask = model.circle_mask(self.i.shape, [5, 5], 4)

    def evaluate(self, q, u, i=None, **kw):
        return model.eye_result(self.i if i is None else i,
                                np.full_like(self.i, q), np.full_like(self.i, u),
                                self.mask, **kw)

    def test_unpolarized_zero_contrast(self):
        result = self.evaluate(0, 0)
        self.assertEqual(result["prc"], 0)
        self.assertIsNone(result["max_axis_deg_mod_180"])

    def test_uniform_sky_analytic_ratio(self):
        result = self.evaluate(.3, .4)
        r, t = 10**(-.4), 10**(-.04)
        b, a = 1-(r+t)/2, (t-r)/2
        self.assertAlmostEqual(result["prc"], math.log((b+a*.5)/(b-a*.5)))
        self.assertAlmostEqual(result["response_max_mean"], b+a*.5)

    def test_analytic_max_matches_independent_pixelwise_heading_sweep(self):
        rng=np.random.default_rng(20261112)
        i=rng.uniform(1,4,(11,11)); dop=rng.uniform(0,.8,(11,11))
        aop=rng.uniform(-math.pi/2,math.pi/2,(11,11))
        q,u=i*dop*np.cos(2*aop),i*dop*np.sin(2*aop)
        result=model.eye_result(i,q,u,self.mask)
        r,t=10**(-.4),10**(-.04)
        logs=[]
        for theta in np.deg2rad(np.arange(0,180,.1)):
            # Independent per-direction DoP/AoP absorption then sum per eye.
            projected=i*dop*np.cos(2*(aop-theta))
            first=i-.5*(t*(i+projected)+r*(i-projected))
            second=i-.5*(t*(i-projected)+r*(i+projected))
            logs.append(math.log(first[self.mask].sum()/second[self.mask].sum()))
        self.assertLessEqual(max(logs),result["prc"]+1e-12)
        self.assertAlmostEqual(max(logs),result["prc"],places=6)

    def test_matches_literal_recovered_R_fixed_axis(self):
        r, t = 10**(-.4), 10**(-.04)
        p0, p1 = 1., .5
        perp = p0-.5*(t*(p0-p1)+r*(p0+p1))
        par = p0-.5*(t*(p0+p1)+r*(p0-p1))
        self.assertAlmostEqual(self.evaluate(p1, 0)["prc"], abs(math.log(par/perp)))

    def test_rotating_global_axis_does_not_change_maximum(self):
        baseline = self.evaluate(.3, .4)
        phi = math.radians(37)
        q, u = .3*math.cos(2*phi)-.4*math.sin(2*phi), .3*math.sin(2*phi)+.4*math.cos(2*phi)
        rotated = self.evaluate(q, u)
        self.assertAlmostEqual(baseline["prc"], rotated["prc"])
        delta = (rotated["max_axis_deg_mod_180"]-baseline["max_axis_deg_mod_180"]) % 180
        self.assertAlmostEqual(delta, 37)

    def test_integrate_before_maximize_opposed_patches_cancel(self):
        mask = np.ones((1, 2), dtype=bool)
        result = model.eye_result(np.ones((1,2)), np.array([[.8,-.8]]),
                                   np.zeros((1,2)), mask)
        self.assertEqual(result["prc"], 0)
        self.assertGreater(self.evaluate(.8, 0)["prc"], 0)

    def test_unpolarized_addition_decreases_contrast(self):
        original = self.evaluate(.3, .4)
        brighter = self.evaluate(.3, .4, i=self.i*2)
        self.assertLess(brighter["prc"], original["prc"])
        self.assertAlmostEqual(original["response_max_sum"]-original["response_min_sum"],
                               brighter["response_max_sum"]-brighter["response_min_sum"])

    def test_intensity_scaling_cancels_and_sums_equal_means(self):
        baseline = self.evaluate(.3, .4)
        scaled = self.evaluate(3, 4, i=self.i*10)
        self.assertAlmostEqual(baseline["prc"], scaled["prc"])
        self.assertAlmostEqual(baseline["prc"], math.log(
            baseline["response_max_mean"]/baseline["response_min_mean"]))

    def test_normalized_direction_is_explicit_distinct_method(self):
        i = np.array([[1., 10.]])
        q = np.array([[.8, 0.]])
        u, mask = np.zeros((1,2)), np.ones((1,2), bool)
        physical = model.eye_result(i, q, u, mask)
        normalized = model.eye_result(i, q, u, mask, response_mode="normalized_direction")
        self.assertGreater(normalized["prc"], physical["prc"])

    def test_invalid_pixels_not_silently_zeroed(self):
        with self.assertRaises(ValueError):
            self.evaluate(2, 0)
        q = np.zeros_like(self.i)
        q[5,5] = 2
        result = model.eye_result(self.i, q, np.zeros_like(q), self.mask,
                                  invalid_policy="exclude")
        self.assertEqual(result["excluded_pixels"], 1)

    def test_bad_geometry_and_absorbances_rejected(self):
        for center, radius in (([5,5], -1), ([0,0], 5), ([5,5], float("nan"))):
            with self.assertRaises(ValueError):
                model.circle_mask((11,11), center, radius)
        with self.assertRaises(ValueError):
            self.evaluate(0, 0, a_parallel=-.1)

    def test_registration_clockwise_and_unsupported_border(self):
        image = np.zeros((9,9)); image[2,4] = 1
        rotated = model.align_image(image, 90, [0,0])
        self.assertAlmostEqual(rotated[4,6], 1)
        shifted = model.align_image(image, 0, [1,0])
        self.assertTrue(np.isnan(shifted[:,0]).all())
        self.assertAlmostEqual(shifted[2,5], 1)

    def test_legacy_pair_normalized_display_rules_preserved(self):
        z, f, n, o = [np.array([[v]]) for v in (.8, .7, .2, .3)]
        i,q,u,qa = model.combine_replicates([[z,f,n,o]]*3, "legacy_dop_aop")
        self.assertAlmostEqual(float(q[0,0]), .6)
        self.assertAlmostEqual(float(u[0,0]), .4)
        z,n = np.ones((1,1)), np.zeros((1,1))
        i,q,u,qa = model.combine_replicates([[z,f,n,o]], "legacy_dop_aop")
        self.assertEqual(float(q[0,0]), 0)
        self.assertEqual(float(u[0,0]), 0)

    def test_legacy_aop_and_linear_stokes_not_silently_interchanged(self):
        # Two polarization angles near opposite ends of the AoP branch cut.
        def quartet(angle):
            q,u = .5*math.cos(2*angle), .5*math.sin(2*angle)
            return [np.full((1,1),v) for v in ((1+q)/2,(1+u)/2,(1-q)/2,(1-u)/2)]
        pairs = [quartet(math.radians(89)),quartet(math.radians(-89))]
        legacy = model.combine_replicates(pairs, "legacy_dop_aop")
        linear = model.combine_replicates(pairs, "linear_stokes")
        self.assertGreater(float(legacy[1][0,0]),0)
        self.assertLess(float(linear[1][0,0]),0)


class FileWorkflowTests(unittest.TestCase):
    def test_discovery_partial_duplicate_and_other_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)/"raw"; root.mkdir()
            for n in range(1,13):
                (root/f"{n}.tif").touch()
            partial=root/"partial"; partial.mkdir(); (partial/"1.tif").touch()
            duplicate=root/"duplicate"; duplicate.mkdir()
            (duplicate/"1.tif").touch(); (duplicate/"01.TIFF").touch()
            other=root/"other"; other.mkdir(); (other/"sky.NEF").touch()
            with contextlib.redirect_stdout(io.StringIO()):
                report=model.discover([root], Path(directory)/"inventory")
            self.assertEqual(sum(c["complete_1_to_12"] for c in report["candidates"]),1)
            self.assertEqual(len(report["other_image_folders"]),1)
            candidate=json.loads((Path(directory)/"inventory/candidate_0000.json").read_text())
            self.assertFalse(candidate["geometry_reviewed"])
            self.assertIsNone(candidate["fov"]["radius_60_px"])

    def make_config(self, path):
        return {"schema_version":1,"geometry_reviewed":True,"registration_reviewed":True,
            "site_id":None,"input":{"kind":"stokes_npz","path":str(path)},
            "fov":{"center_xy_px":[5,5],"radius_60_px":4,"angular_radius_deg":60},
            "absorbance_parallel":.4,"absorbance_perpendicular":.04,
            "response_mode":"stokes_intensity","invalid_policy":"error"}

    def test_npz_end_to_end_raw_preservation_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/"sky.npz"
            np.savez(source,I=np.ones((11,11)),Q=np.full((11,11),.3),U=np.full((11,11),.4))
            before=model.digest(source)
            config=self.make_config(source); manifest=root/"config.json"
            manifest.write_text(json.dumps(config))
            out=root/"result"
            with contextlib.redirect_stdout(io.StringIO()):
                result=model.run(manifest,out)
            self.assertGreater(result["prc"],0)
            self.assertEqual(before,model.digest(source))
            self.assertTrue((out/"qc_preview.png").is_file())
            with self.assertRaises(FileExistsError):
                model.run(manifest,out)

    def test_unreviewed_and_unknown_radius_fail_before_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); manifest=root/"config.json"
            config=self.make_config(root/"missing.npz")
            config["geometry_reviewed"]=False
            manifest.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError,"Review geometry"):
                model.run(manifest,root/"out")
            config["fov"]["radius_60_px"]=None
            manifest.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError,"measured"):
                model.run(manifest,root/"out",preview_only=True)
            self.assertFalse((root/"out").exists())

    def test_dop_aop_degrees_and_radians_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            outputs=[]
            for unit,angle in (("degrees",30.),("radians",math.pi/6)):
                path=root/f"{unit}.npz"
                np.savez(path,DoP=np.full((11,11),.5),AoP=np.full((11,11),angle))
                config=self.make_config(path)
                config["input"]={"kind":"dop_aop_npz","path":str(path),"aop_unit":unit}
                config["response_mode"]="normalized_direction"
                outputs.append(model.load_sky(config,root)[:3])
            for left,right in zip(*outputs):
                np.testing.assert_allclose(left,right)

    def test_dop_maps_require_intensity_or_explicit_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);path=root/"sky.npz"
            np.savez(path,DoP=np.ones((11,11))*.2,AoP=np.zeros((11,11)))
            config=self.make_config(path)
            config["input"]={"kind":"dop_aop_npz","path":str(path),"aop_unit":"radians"}
            with self.assertRaisesRegex(ValueError,"must include I"):
                model.load_sky(config,root)

    @unittest.skipUnless(importlib.util.find_spec("tifffile"),"optional TIFF dependency")
    def test_tiff_triplets_keep_uint16_pixels_and_match_known_sky(self):
        import tifffile
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); raw=root/"raw"; raw.mkdir()
            # Integer analyzer intensities yield I=1000,Q=300,U=400.
            values=[650,700,350,300]
            for group in model.GROUPS:
                for n,value in zip(group,values):
                    tifffile.imwrite(raw/f"{n}.tif",np.full((11,11,3),value,dtype=np.uint16))
            config=self.make_config(raw)
            config["input"]={"kind":"numbered_tiffs","directory":str(raw),"channel":2,
                "groups":[[f"{n}.tif" for n in g] for g in model.GROUPS],
                "rotation_deg_clockwise":[0,0,0,0],"shift_xy_px":[[0,0]]*4,
                "stokes_method":"linear_stokes"}
            hashes={p.name:model.digest(p) for p in raw.iterdir()}
            i,q,u,records,qa=model.load_sky(config,root)
            self.assertTrue(np.all(i==1000));self.assertTrue(np.all(q==300));self.assertTrue(np.all(u==400))
            self.assertEqual(len(records),12)
            self.assertEqual(records[0]["dtype"],"uint16")
            self.assertEqual(hashes,{p.name:model.digest(p) for p in raw.iterdir()})

    @unittest.skipUnless(importlib.util.find_spec("tifffile"),"optional TIFF dependency")
    def test_single_explicitly_ordered_tiff_replicate(self):
        import tifffile
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);raw=root/"raw";raw.mkdir()
            names=[f"{n}.tif" for n in range(1,5)]
            for name,value in zip(names,[650,700,350,300]):
                tifffile.imwrite(raw/name,np.full((11,11,3),value,dtype=np.uint16))
            config=self.make_config(raw)
            config["input"]={"kind":"tiff_groups","directory":str(raw),"channel":2,
                "groups":[names],"rotation_deg_clockwise":[0,0,0,0],
                "shift_xy_px":[[0,0]]*4,"stokes_method":"linear_stokes"}
            i,q,u,records,qa=model.load_sky(config,root)
            np.testing.assert_array_equal(i,np.full((11,11),1000))
            np.testing.assert_array_equal(q,np.full((11,11),300))
            np.testing.assert_array_equal(u,np.full((11,11),400))
            self.assertEqual(len(qa),1);self.assertEqual(len(records),4)
            # These synthetic filenames carry a TEST-defined order; no field order inferred.
            config["input"]["groups"]=[names,names,names]
            with self.assertRaisesRegex(ValueError,"do not duplicate"):
                model.load_sky(config,root)

    def test_discovery_marks_four_images_candidate_not_three_replicates(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);raw=root/"raw";raw.mkdir()
            for n in range(1,5):
                (raw/f"{n}.tif").touch()
            with contextlib.redirect_stdout(io.StringIO()):
                result=model.discover([raw],root/"inventory")
            entry=result["candidates"][0]
            self.assertTrue(entry["single_1_to_4_candidate"])
            self.assertFalse(entry["complete_1_to_12"])
            config=json.loads((root/"inventory/candidate_0000.json").read_text())
            self.assertEqual(config["input"]["groups"],[["1.tif","2.tif","3.tif","4.tif"]])
            self.assertFalse(config["registration_reviewed"])

    @unittest.skipUnless(importlib.util.find_spec("tifffile"),"optional TIFF dependency")
    def test_header_inspection_and_unaligned_quicklook(self):
        import tifffile
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);raw=root/"raw";raw.mkdir()
            names=[f"DSC_{n}.tif" for n in range(5080,5084)]
            for name in names:
                tifffile.imwrite(raw/name,np.full((11,11,3),100,dtype=np.uint16))
            inventory=root/"inventory.json"
            inventory.write_text(json.dumps({"candidates":[],"other_image_folders":[
                {"folder":str(raw),"first_20_filenames":names}]}))
            before={p.name:model.digest(p) for p in raw.iterdir()}
            with contextlib.redirect_stdout(io.StringIO()):
                report=model.inspect_inventory(inventory,root/"inspection",raw)
            self.assertEqual(report["tiff_headers_read"],4)
            self.assertEqual(report["errors"],[])
            self.assertEqual(report["records"][0]["page0_dtype"],"uint16")
            self.assertTrue((root/"inspection/input_quicklook.png").is_file())
            self.assertEqual(before,{p.name:model.digest(p) for p in raw.iterdir()})

    @unittest.skipUnless(importlib.util.find_spec("tifffile"),"optional TIFF dependency")
    def test_header_inspection_reports_unreadable_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);raw=root/"raw";raw.mkdir();(raw/"1.tif").write_text("bad file")
            inventory=root/"inventory.json"
            inventory.write_text(json.dumps({"candidates":[{"folder":str(raw),"files":{"1":["1.tif"]}}]}))
            with contextlib.redirect_stdout(io.StringIO()):
                report=model.inspect_inventory(inventory,root/"inspection")
            self.assertEqual(report["tiff_headers_read"],0)
            self.assertEqual(len(report["errors"]),1)


if __name__=="__main__":
    unittest.main()
