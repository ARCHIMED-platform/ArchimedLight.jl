import java.util.ArrayList;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3f;
import fr.amap.archimed.application.process.voxel.Ray;
import fr.amap.archimed.application.process.voxel.RayPath;
import fr.amap.archimed.application.process.voxel.VoxelInfos;
import fr.amap.archimed.application.process.voxel.VoxelSpace;
import fr.amap.archimed.application.process.voxel.VoxelSpaceInfos;

public class VoxelBaseline {
    private static VoxelSpace grid(int nx, int ny, int nz, double sx, double sy, double sz, float pad) {
        VoxelSpace vs = new VoxelSpace(new VoxelSpaceInfos(
            new Point3d(0, 0, 0), new Point3d(sx, sy, sz), new Point3i(nx, ny, nz)));
        vs.voxels = new ArrayList<>();
        for (int n = 0; n < nx * ny * nz; n++) {
            VoxelInfos voxel = new VoxelInfos(n);
            voxel.PAD = pad;
            vs.voxels.add(voxel);
        }
        return vs;
    }

    private static void path(String name, VoxelSpace vs, Point3d origin, Vector3f direction) {
        direction.normalize();
        Ray ray = new Ray(vs, origin);
        ray.crossVoxels(direction);
        System.out.println("PATH\t" + name + "\t" + direction.x + "\t" + direction.y + "\t" + direction.z);
        for (RayPath segment : ray.getCrossedPaths()) {
            int index = segment.crossedVoxel == null ? -1 : segment.crossedVoxel.listIndex;
            System.out.printf(java.util.Locale.ROOT, "SEGMENT\t%s\t%d\t%.17g%n",
                name, index, segment.pathLength);
        }
    }

    private static void response(String name, boolean toric, Vector3f direction) {
        VoxelSpace vs = grid(3, 2, 3, 3, 2, 3, 0.4f);
        vs.setRequestedRayCount(4);
        vs.setDirectionCount(1);
        direction.normalize();
        vs.computeInterceptionCoefForDir(0, direction, 0.5f, toric);
        vs.computeInterceptedIntensity(0, 2.0f);
        double total = 0;
        for (VoxelInfos voxel : vs.voxels) {
            total += voxel.intercepted;
            System.out.printf(java.util.Locale.ROOT, "VOXEL\t%s\t%d\t%.17g\t%.17g%n",
                name, voxel.listIndex, voxel.interceptionCoefs[0], voxel.intercepted);
        }
        System.out.printf(java.util.Locale.ROOT, "TOTAL\t%s\t%.17g\t%d%n",
            name, total, vs.getEffectiveRayCount());
    }

    private static void profileResponse(
        String name, VoxelSpace vs, int rayCount, boolean toric, Vector3f direction
    ) {
        vs.setRequestedRayCount(rayCount);
        vs.setDirectionCount(1);
        direction.normalize();
        vs.computeInterceptionCoefForDir(0, direction, 0.5f, toric);
        vs.computeInterceptedIntensity(0, 2.0f);

        int nz = vs.getVoxelSpaceInfos().getSplit().z;
        double[] coefficients = new double[nz];
        double[] intercepted = new double[nz];
        double total = 0;
        for (VoxelInfos voxel : vs.voxels) {
            int k = voxel.listIndex % nz;
            coefficients[k] += voxel.interceptionCoefs[0];
            intercepted[k] += voxel.intercepted;
            total += voxel.intercepted;
        }
        for (int k = 0; k < nz; k++) {
            System.out.printf(java.util.Locale.ROOT, "PROFILE\t%s\t%d\t%.17g\t%.17g%n",
                name, k, coefficients[k], intercepted[k]);
        }
        System.out.printf(java.util.Locale.ROOT, "TOTAL\t%s\t%.17g\t%d%n",
            name, total, vs.getEffectiveRayCount());
    }

    public static void main(String[] args) {
        VoxelSpace pathGrid = grid(2, 2, 3, 2, 2, 3, 0.4f);
        path("vertical", pathGrid, new Point3d(0.5, 0.5, 3), new Vector3f(0, 0, -1));
        path("diagonal", pathGrid, new Point3d(0.25, 0.5, 3), new Vector3f(0.5f, 0, -1));
        path("corner", pathGrid, new Point3d(0.5, 0.5, 3), new Vector3f(1, 1, -1));

        response("vertical_toric", true, new Vector3f(0, 0, -1));
        response("diagonal_toric", true, new Vector3f(0.5f, 0.2f, -1));
        response("diagonal_java_nontoric", false, new Vector3f(0.5f, 0.2f, -1));

        VoxelSpace ground = grid(5, 5, 2, 5, 5, 2, 0.0f);
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                ground.voxels.get(i * 5 * 2 + j * 2).PAD = 0.5f;
        profileResponse("ground_toric", ground, 4, true, new Vector3f(0.5f, 0.2f, -1));

        VoxelSpace pillar = grid(5, 5, 4, 5, 5, 4, 0.0f);
        for (int k = 0; k < 4; k++)
            pillar.voxels.get(2 * 5 * 4 + 2 * 4 + k).PAD = 0.5f;
        profileResponse(
            "pillar_java_nontoric", pillar, 9, false, new Vector3f(0.05f, 0.0f, -1)
        );
    }
}
