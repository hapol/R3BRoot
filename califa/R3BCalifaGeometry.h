#ifndef R3BCALIFAGEOMETRY_H
#define R3BCALIFAGEOMETRY_H

#include <TObject.h>

class TVector3;
class TGeoNavigator;

/**
 * Geometrical queries to CALIFA
 */
class R3BCalifaGeometry : public TObject
{
  public:
    /** Default constructor */
    R3BCalifaGeometry();

    /** Standard constructor.
     *@param version geonetry version
     */
    R3BCalifaGeometry(int version);

    /** Destructor **/
    ~R3BCalifaGeometry();

    /**
     * Gets position in polar coordinates of crystal with given ID.
     *
     * @param iD crystal ID (depending on geometry version)
     * @param polar [out] Will be filled with polar angle (radians) of crystal center
     * @param azimuthal [out] Will be filled with azimuthal angle (radians) of crystal center
     * @param rho [out] Will be filled with distance (cm) of crystal center to target position (0,0,0)
     */
    Bool_t GetAngles(Int_t iD, Double_t& polar, Double_t& azimuthal, Double_t& rho);

    /**
     * Gets volume path of crystal with given ID.
     *
     * @param iD crystal ID (depending on geometry version)
     * @return Volume path
     */
    const char* GetCrystalVolumePath(int iD);

    /**
     * Gets crystal ID for given volume path.
     *
     * @param volumePath Volume path
     * @return Crystal ID
     */
    int GetCrystalId(const char* volumePath);

    /**
     * Calculate the distance of a given straight track through the active detector volume (crystal(s)). Usefull for
     * iPhos.
     *
     * @param startVertex Start position of track (e.g. (0,0,0) for center of target)
     * @param direction Direction of the track (does not need to be normalized)
     * @param hitPos [out] (optional) If given, this TVector3 will be filled with the position of the first contact
     * between the track and the first crystal along the path.
     * @param numCrystals [in,out]
     * \parblock
     *  (optional) If given together with crystalIDs,
     *  this function will fill numCrystals with the number of crystals along the straight
     *  and the array crystalIds with the IDs of the crystals along the straight.
     *  The array crystalIds MUST be allocated by the caller.
     *  The maximum size of the array has to be given with the content of numCrystals.
     *  The function will stop filling the array after reaching the specified, maximum capacity.
     * \endparblock
     * @param crystalIds [out] (optional, array) See numCrystals.
     * @return Distance (cm) through crystal(s)
     */
    double GetDistanceThroughCrystals(TVector3& startVertex,
                                      TVector3& direction,
                                      TVector3* hitPos = NULL,
                                      int* numCrystals = NULL,
                                      int* crystalIds = NULL);

    /**
     * Returns singleton instance of R3BCalifaGeometry for given geometry version.
     * Warning: Switching the geometry version at run time will lead to undefined beheaviour!
     *
     * @param version Geometry version to use. If in doubt, use 2020.
     * @return Instance of R3BCalifaGeometry
     */
    static R3BCalifaGeometry* Instance(int version);

  private:
    int fGeometryVersion;
    static R3BCalifaGeometry* inst;

    ClassDef(R3BCalifaGeometry, 5);
};

#endif
