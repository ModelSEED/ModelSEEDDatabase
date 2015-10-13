Human metabolic model reconstructed from Recon2.04.mat from VMH.


Methods:

v1.1

Parse reactions for compounds against human metabolome database (hmdb.ca)
Search each term against published literature and open internet (Google Scholar, Google respectively)
Terms that align with plants are eliminated. Terms that align with human are marked human only. Terms that are mycobacterial are marked mycobacterial only. Terms that are microbial but not clearly mycobacterial are eliminated if the microbe or the biochemical process does not involve human conditions (e.g deep sea vents, etc). Otherwise if microbial and unsure if mycobacterial, leave in.

Eliminate: alpha-terpineol, alpha cedrine, beta cucurmene, borneol, camphor, fenchol, ephedrine, "Isopiperitenone" Limonene linelol, menthol maackiain methone sophorol thujan Carveol Pinene vetispiradiene santalene caryophyllene  camphane sabinone Sabinene columbianetin Dihydrocarvone"  Isozizaene" Germacrene  "Pulegone" "Acenaphthoquinone" "Nigerose" Sesquiphellandrene Geosmin Eudesmol"  Selinene" "dihydroxicinnamic acid" Muurola-3,5-diene - Bing
 "Humulene" Guaiene "Avermitilol" - "Longifolene" "Viridiflorene" "Zingiberene" Bisabolene" Capreomycidine" Bergamotene" "Eriodictyol Dihydrozeatin zeatin Daturine[0] Brassica Solasodine[0] Stigmasterol[0] Vitexin[0] coniferol yeast-specific  yersiniose Taxifolin[0] Eriodictyol[0] Luteoforol[0] Aureusidin[0] Bracteatin[0] Naringenin[0] achilleol  Amyrin[0] Avenacin Baccharis Baruol[0] Camelliol Friedelin[0]  Vinblastine[0] Vindoline[0] + (1) Catharanthine[0  Horhammericine[0]  Tabersonine[0] Lochnerinine  Minovincinine[0] Neoxanthin[0]  Violaxanthin[0] Xanthoxin[0] Abscisic Butein[0] Ampelopsin[0] Pinocembrin[0] Liquiritigenin[0]  Flavanone coniferonic Afzelechin[0] Raucaffricine[0] Humulene[0] Nerolidol[0] Euphol Azadirachtin[0] Cucurbitacin Limonoate[0]

Glucogallin[0] Fenchyl Sabinene Zeatin[0] Quercetin[0] Cadinene[0] Isopiperitenol[0] Perillyl Prunin[0] Hydroxychalcone[0] Camphor[0] Borneol[0] Gibberellin Apigenin[0] Isosalipurpol[0] Resveratrol[0] Piceatannol[0] Protocatechuate[0]  zeatin quercetin kanamycin Caffeoyl scyllo-Inosamine[0] Orcinol[0] Vanillin[0] linatine  Cycloartenol[0] Desulfoglucotropeolin[0] Glyceollin[0] Sequoyitol[0] Loganin[0] secloganin (override secloganin+dopamine...eliminate secloganin) Salicin[0] Pectin[0] Inulin[0] Stachyose[0] Scopoletin[0] Geranial[0] Syringin[0]
 Vicianin[0] Stipitatate[0] Phytate[0] Dipicolinate[0]  Stachydrine[0] Hypoglycine betanidin Zeacarotene[0] Cryptoxanthin Violaxanthin[0] Xanthoxin[0] Lotaustralin[0] phytol chlorophyll

Uneliminated: Catechol (found papers on human and mycobaterium metabolism)

Human:  UDP-N-acetylglucosamine[0] Renin Carnosine[0] (Bhatt et al suggests smegmatis can use carnosine) Prostaglandin "Factor" Leukotriene Complement  Dopaquinone[0] creatinine carnitine Sphingosine Psychosine[0] GABA[0] Dopa Calcidiol[0] adrenodoxin[0] Dermatan[0] Hyaluronate[0]

Myco:  UDP-N-acetylglucosamine[0] Neurosporene[0] Spheroidene[0] teichoic acid (note trehalose 2 sulfate reaction found (rxn10671)

? alpha terpineol, Glucocerebroside (may be in both?) Melatonin[0] (found in animals/plants/fungi/bacteria)

v1.0
XRemoved all reactions with plastids, thylakoids in name
XRemoved Gibberellins.
XParse substrates against MetaCyc human 



BiomassCompounds

#MATLAB commands used to parse Recon2 for objective function. Find reaction matching pattern biomass, extract S matrix components, then extract names for conversion to ModelSEED compounds.

>> strmatch('biomass',modelR204.rxns)

ans =

        6799

>> modelR204.rxns(6799)

ans = 

    'biomass_reaction'

modelR204.S(:,6799)   #This command checks the column of the S matrix corresponding to the biomass reaction

  (20,1)     -20.6508
  (37,1)     -20.7045
  (39,1)      20.6508
  (40,1)      20.6508
  (41,1)      20.6508
 (126,1)      -0.3859
 (149,1)      -0.3526
 (437,1)      -0.0361
 (571,1)      -0.2794
 (572,1)      -0.5056
 (574,1)      -0.0466
 (577,1)      -0.3260
 (580,1)      -0.5389
 (582,1)      -0.3925
 (585,1)      -0.3127
 (648,1)      -0.5921
 (679,1)      -0.3593
 (824,1)      -0.1530
 (923,1)      -0.0233
 (925,1)      -0.0390
 (929,1)      -0.1545
 (930,1)      -0.0554
 (951,1)      -0.0204
 (983,1)      -0.0029
 (984,1)      -0.0117
(1070,1)      -0.0534
(1101,1)      -0.0099
(1109,1)      -0.0094
(1110,1)      -0.0132
(1240,1)      -0.0131
(1860,1)      -0.2752
(2112,1)      -0.1264
(2174,1)      -0.1597
(2184,1)      -0.2861
(2222,1)      -0.5455
(2257,1)      -0.0133
(2509,1)      -0.2595
(2558,1)      -0.4125
(2568,1)      -0.0058
(2660,1)      -0.0175
(2750,1)      -0.3526

modelR204.metNames(20)
modelR204.metNames(37)
modelR204.metNames(39)
modelR204.metNames(40)
modelR204.metNames(41)
modelR204.metNames(126)
modelR204.metNames(149)
modelR204.metNames(437)
modelR204.metNames(571)
modelR204.metNames(572)
modelR204.metNames(574)
modelR204.metNames(577)
modelR204.metNames(580)
modelR204.metNames(582)
modelR204.metNames(585)
modelR204.metNames(648)
modelR204.metNames(679)
modelR204.metNames(824)
modelR204.metNames(923)
modelR204.metNames(925)
modelR204.metNames(929)
modelR204.metNames(930)
modelR204.metNames(951)
modelR204.metNames(983:984)
modelR204.metNames(1070)
modelR204.metNames(1101)
modelR204.metNames(1109:1110)
modelR204.metNames(1240)
modelR204.metNames(1860)
modelR204.metNames(2112)
modelR204.metNames(2174)
modelR204.metNames(2184)
modelR204.metNames(2222)
modelR204.metNames(2257)
modelR204.metNames(2509)
modelR204.metNames(2558)
modelR204.metNames(2568)
modelR204.metNames(2660)
modelR204.metNames(2750)
    'Water'
    'ATP'
    'ADP'
    'proton'
    'hydrogenphosphate'
    'L-glutamate(1-)'
    'L-aspartate(1-)'
    'GTP'
    'L-asparagine'
    'L-alanine'
    'L-cysteine'
    'L-glutamine'
    'Glycine'
    'L-serine'
    'L-threonine'
    'L-lysinium(1+)'
    'L-argininium(1+)'
    'L-methionine'
    '1-phosphatidyl-1D-myo-inositol(1-)'
    'CTP'
    'Phosphatidylcholine'
    'phosphatidylethanolamine'
    'cholesterol'
    'phosphatidylglycerol(1-)'
    'cardiolipin'
    'UTP'
    'dGTP'
    'dCTP'
    'dATP'
    'dTTP'
    'D-Glucose 6-phosphate'
    'L-histidine'
    'L-tyrosine'
    'L-isoleucine'
    'L-leucine'
    'L-tryptophan'
    'L-phenylalanine'
    'L-proline'
    'phosphatidylserine'
    'sphingomyelin betaine'
    'L-valine'