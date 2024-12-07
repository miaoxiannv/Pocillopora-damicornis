"""
读取mapping文件，找到source列中为overlapping的对应的protein_id,输出结果为列表基因的
protein_id	source
XP_022779499_1	marker1_only
XP_022779618_1	marker1_only
XP_022779926_1	marker1_only
XP_022780118_1	marker1_only
XP_022781572_1	marker1_only
"""

def read_mapping_file(mapping_files):
    """读取多个mapping文件，返回overlapping的protein_id和对应文件名的字典"""
    overlapping_proteins_dict = {}  # 改用字典存储蛋白质ID和来源文件的关系
    
    # 处理单个文件路径或文件路径列表
    if isinstance(mapping_files, str):
        mapping_files = [mapping_files]
        
    for mapping_file in mapping_files:
        file_name = mapping_file.split('\\')[-1].replace('.tsv', '')  # 获取文件名（不含路径和扩展名）
        with open(mapping_file, 'r') as f:
            # 跳过表头
            next(f)
            for line in f:
                protein_id, source = line.strip().split('\t')
                if source == 'overlapping':
                    if protein_id not in overlapping_proteins_dict:
                        overlapping_proteins_dict[protein_id] = []
                    overlapping_proteins_dict[protein_id].append(file_name)
    
    return overlapping_proteins_dict


"""

读取两个物种orthofinder结果，寻找出1：1的结果，高保守型基因
建立
Orthogroup	FL	SP-cell
OG0000000		"Spis10168_1, Spis10297_1, Spis10529_1, Spis10717_1, Spis12045_1, Spis13206_1, Spis13628_1, Spis13647_1, Spis13688_1, Spis14197_1, Spis14215_1, Spis14563_1, Spis14568_1, Spis14757_1, Spis14770_1, Spis15806_1, Spis16318_1, Spis16834_1, Spis16968_1, Spis17533_1, Spis17772_1, Spis18071_1, Spis19035_1, Spis19087_1, Spis19421_1, Spis19605_1, Spis1977_1, Spis1979_1, Spis19830_1, Spis19831_1, Spis1985_1, Spis1987_1, Spis20355_1, Spis20624_1, Spis2090_1, Spis20975_1, Spis21017_1, Spis21018_1, Spis21034_1, Spis21625_1, Spis21626_1, Spis21795_1, Spis21798_1, Spis21942_1, Spis22013_1, Spis22682_1, Spis22811_1, Spis22988_1, Spis23022_1, Spis23034_1, Spis23089_1, Spis23121_1, Spis23122_1, Spis23285_1, Spis23586_1, Spis23660_1, Spis23662_1, Spis23667_1, Spis23891_1, Spis24059_1, Spis24074_1, Spis24210_1, Spis24361_1, Spis24490_1, Spis24687_1, Spis24886_1, Spis25265_1, Spis25270_1, Spis252_1, Spis25554_1, Spis25606_1, Spis25648_1, Spis25725_1, Spis2594_1, Spis333_1, Spis3444_1, Spis3552_1, Spis3557_1, Spis4366_1, Spis5247_1, Spis6072_1, Spis6073_1, Spis6600_1, Spis6728_1, Spis7509_1, Spis8267_1, Spis8307_1, Spis8469_1, Spis8655_1, Spis8947_1, Spis9441_1, Spis9564_1, Spis9990_1, Spis9994_1, Spis_XP_022777559_1, Spis_XP_022777768_1, Spis_XP_022777857_1, Spis_XP_022777858_1, Spis_XP_022778085_1, Spis_XP_022778196_1, Spis_XP_022778259_1, Spis_XP_022778282_1, Spis_XP_022778283_1, Spis_XP_022778379_1, Spis_XP_022778796_1, Spis_XP_022778797_1, Spis_XP_022779006_1, Spis_XP_022779007_1, Spis_XP_022779457_1, Spis_XP_022779469_1, Spis_XP_022780505_1, Spis_XP_022781742_1, Spis_XP_022781755_1, Spis_XP_022782764_1, Spis_XP_022783263_1, Spis_XP_022783264_1, Spis_XP_022783298_1, Spis_XP_022783913_1, Spis_XP_022784027_1, Spis_XP_022784343_1, Spis_XP_022784574_1, Spis_XP_022784631_1, Spis_XP_022784654_1, Spis_XP_022784666_1, Spis_XP_022785230_1, Spis_XP_022785231_1, Spis_XP_022786152_1, Spis_XP_022786310_1, Spis_XP_022787478_1, Spis_XP_022787688_1, Spis_XP_022787689_1, Spis_XP_022787690_1, Spis_XP_022788132_1, Spis_XP_022788133_1, Spis_XP_022788533_1, Spis_XP_022789237_1, Spis_XP_022790232_1, Spis_XP_022790233_1, Spis_XP_022790440_1, Spis_XP_022790960_1, Spis_XP_022790962_1, Spis_XP_022790963_1, Spis_XP_022790964_1, Spis_XP_022790966_1, Spis_XP_022791300_1, Spis_XP_022791541_1, Spis_XP_022791786_1, Spis_XP_022791844_1, Spis_XP_022791987_1, Spis_XP_022794233_1, Spis_XP_022794387_1, Spis_XP_022795007_1, Spis_XP_022795008_1, Spis_XP_022795598_1, Spis_XP_022795745_1, Spis_XP_022796094_1, Spis_XP_022796432_1, Spis_XP_022796861_1, Spis_XP_022797506_1, Spis_XP_022797507_1, Spis_XP_022797520_1, Spis_XP_022797621_1, Spis_XP_022797977_1, Spis_XP_022798022_1, Spis_XP_022798305_1, Spis_XP_022799297_1, Spis_XP_022799861_1, Spis_XP_022799869_1, Spis_XP_022799870_1, Spis_XP_022799900_1, Spis_XP_022800628_1, Spis_XP_022801004_1, Spis_XP_022801248_1, Spis_XP_022801341_1, Spis_XP_022801565_1, Spis_XP_022801986_1, Spis_XP_022802575_1, Spis_XP_022802779_1, Spis_XP_022802837_1, Spis_XP_022803018_1, Spis_XP_022804234_1, Spis_XP_022804327_1, Spis_XP_022804737_1, Spis_XP_022804765_1, Spis_XP_022804766_1, Spis_XP_022804786_1, Spis_XP_022805148_1, Spis_XP_022805184_1, Spis_XP_022805284_1, Spis_XP_022805285_1, Spis_XP_022805463_1, Spis_XP_022805940_1, Spis_XP_022805942_1, Spis_XP_022805950_1, Spis_XP_022805992_1, Spis_XP_022806184_1, Spis_XP_022806462_1, Spis_XP_022806477_1, Spis_XP_022806557_1, Spis_XP_022807191_1, Spis_XP_022807271_1, Spis_XP_022807272_1, Spis_XP_022807382_1, Spis_XP_022807383_1, Spis_XP_022807985_1, Spis_XP_022808364_1, Spis_XP_022808550_1, Spis_XP_022808694_1, Spis_XP_022808894_1, Spis_XP_022808943_1, Spis_XP_022808944_1, Spis_XP_022809212_1, Spis_XP_022809213_1, Spis_XP_022809281_1, Spis_XP_022809470_1, Spis_XP_022809517_1, Spis_XP_022809686_1, Spis_XP_022809721_1, Spis_XP_022809970_1, Spis_XP_022810113_1, Spis_XP_022810133_1, Spis_XP_022810538_1"
OG0000001		"Spis10012_1, Spis10075_1, Spis10078_1, Spis10187_1, Spis10190_1, Spis10191_1, Spis10319_1, Spis10451_1, Spis10465_1, Spis10466_1, Spis10487_1, Spis10619_1, Spis10627_1, Spis10634_1, Spis11004_1, Spis11233_1, Spis11664_1, Spis12263_1, Spis1230_1, Spis12318_1, Spis12406_1, Spis12466_1, Spis12469_1, Spis12638_1, Spis12885_1, Spis13077_1, Spis13141_1, Spis13179_1, Spis13326_1, Spis13348_1, Spis13433_1, Spis13715_1, Spis14144_1, Spis14308_1, Spis14432_1, Spis14517_1, Spis14643_1, Spis14765_1, Spis14875_1, Spis14902_1, Spis14948_1, Spis14971_1, Spis14976_1, Spis15390_1, Spis15455_1, Spis15722_1, Spis15783_1, Spis15943_1, Spis16535_1, Spis16691_1, Spis16837_1, Spis16923_1, Spis17252_1, Spis17259_1, Spis17611_1, Spis17790_1, Spis17812_1, Spis18206_1, Spis18282_1, Spis18343_1, Spis18680_1, Spis18746_1, Spis19076_1, Spis19134_1, Spis19340_1, Spis19388_1, Spis19442_1, Spis19848_1, Spis20067_1, Spis2007_1, Spis20131_1, Spis20242_1, Spis20280_1, Spis202_1, Spis20326_1, Spis20360_1, Spis20548_1, Spis20651_1, Spis20998_1, Spis21026_1, Spis2113_1, Spis21458_1, Spis21566_1, Spis21859_1, Spis21963_1, Spis22253_1, Spis22315_1, Spis22389_1, Spis22520_1, Spis22623_1, Spis22822_1, Spis22824_1, Spis22934_1, Spis23008_1, Spis23015_1, Spis23032_1, Spis23037_1, Spis23162_1, Spis2319_1, Spis23304_1, Spis23308_1, Spis23346_1, Spis23348_1, Spis23463_1, Spis23494_1, Spis23505_1, Spis23600_1, Spis23756_1, Spis23844_1, Spis24049_1, Spis24253_1, Spis24297_1, Spis24438_1, Spis24722_1, Spis24772_1, Spis24912_1, Spis25004_1, Spis25037_1, Spis25515_1, Spis25537_1, Spis25667_1, Spis25700_1, Spis25820_1, Spis2585_1, Spis332_1, Spis3436_1, Spis3891_1, Spis3893_1, Spis396_1, Spis4946_1, Spis6221_1, Spis6301_1, Spis6682_1, Spis687_1, Spis6942_1, Spis7002_1, Spis7016_1, Spis7513_1, Spis7577_1, Spis7786_1, Spis790_1, Spis791_1, Spis8196_1, Spis8284_1, Spis8670_1, Spis8688_1, Spis9357_1, Spis9565_1, Spis9881_1, Spis_XP_022777699_1, Spis_XP_022778374_1, Spis_XP_022779672_1, Spis_XP_022779838_1, Spis_XP_022780456_1, Spis_XP_022781135_1, Spis_XP_022781485_1, Spis_XP_022781959_1, Spis_XP_022782228_1, Spis_XP_022783320_1, Spis_XP_022784874_1, Spis_XP_022785980_1, Spis_XP_022787417_1, Spis_XP_022787862_1, Spis_XP_022788461_1, Spis_XP_022790226_1, Spis_XP_022791346_1, Spis_XP_022791802_1, Spis_XP_022792245_1, Spis_XP_022792420_1, Spis_XP_022793003_1, Spis_XP_022793558_1, Spis_XP_022794170_1, Spis_XP_022794305_1, Spis_XP_022795148_1, Spis_XP_022797595_1, Spis_XP_022800279_1, Spis_XP_022801800_1, Spis_XP_022801966_1, Spis_XP_022802718_1, Spis_XP_022803206_1, Spis_XP_022803272_1, Spis_XP_022804740_1, Spis_XP_022805458_1, Spis_XP_022808140_1, Spis_XP_022808253_1, Spis_XP_022809244_1, Spis_XP_022809752_1, Spis_XP_022810263_1"
OG0000002	transcript_HQ_P2_transcript6075_f7p0_3029_1	"Spis10329_1, Spis1044_1, Spis10471_1, Spis10787_1, Spis12326_1, Spis12413_1, Spis13098_1, Spis13099_1, Spis13325_1, Spis13641_1, Spis14029_1, Spis14234_1, Spis14270_1, Spis14306_1, Spis14410_1, Spis14534_1, Spis16428_1, Spis16462_1, Spis16463_1, Spis16975_1, Spis17165_1, Spis17774_1, Spis18751_1, Spis18818_1, Spis18822_1, Spis18952_1, Spis19096_1, Spis1988_1, Spis20002_1, Spis20060_1, Spis20075_1, Spis20107_1, Spis20364_1, Spis20384_1, Spis20385_1, Spis20727_1, Spis2079_1, Spis21033_1, Spis21117_1, Spis21155_1, Spis21382_1, Spis21444_1, Spis21459_1, Spis2155_1, Spis21725_1, Spis21764_1, Spis21860_1, Spis22061_1, Spis22081_1, Spis22083_1, Spis22176_1, Spis22212_1, Spis22414_1, Spis22635_1, Spis22772_1, Spis22869_1, Spis22908_1, Spis22931_1, Spis22982_1, Spis23006_1, Spis23098_1, Spis23206_1, Spis23294_1, Spis23602_1, Spis23680_1, Spis23871_1, Spis24158_1, Spis24187_1, Spis24309_1, Spis24431_1, Spis24449_1, Spis24706_1, Spis25218_1, Spis25461_1, Spis3235_1, Spis4892_1, Spis4973_1, Spis5148_1, Spis5801_1, Spis5802_1, Spis5947_1, Spis6360_1, Spis6634_1, Spis679_1, Spis7525_1, Spis8549_1, Spis8803_1, Spis928_1, Spis939_1, Spis9442_1, Spis947_1, Spis9532_1, Spis_XP_022777877_1, Spis_XP_022779210_1, Spis_XP_022779241_1, Spis_XP_022779776_1, Spis_XP_022779919_1, Spis_XP_022780500_1, Spis_XP_022781426_1, Spis_XP_022783273_1, Spis_XP_022783300_1, Spis_XP_022783322_1, Spis_XP_022784883_1, Spis_XP_022785142_1, Spis_XP_022786418_1, Spis_XP_022787546_1, Spis_XP_022788128_1, Spis_XP_022788129_1, Spis_XP_022788169_1, Spis_XP_022788174_1, Spis_XP_022790234_1, Spis_XP_022790318_1, Spis_XP_022791503_1, Spis_XP_022791845_1, Spis_XP_022791847_1, Spis_XP_022791988_1, Spis_XP_022792295_1, Spis_XP_022792296_1, Spis_XP_022792658_1, Spis_XP_022792659_1, Spis_XP_022793396_1, Spis_XP_022793931_1, Spis_XP_022795609_1, Spis_XP_022796124_1, Spis_XP_022796419_1, Spis_XP_022796430_1, Spis_XP_022796469_1, Spis_XP_022796583_1, Spis_XP_022797032_1, Spis_XP_022797033_1, Spis_XP_022797448_1, Spis_XP_022797518_1, Spis_XP_022797622_1, Spis_XP_022797858_1, Spis_XP_022797971_1, Spis_XP_022800415_1, Spis_XP_022800430_1, Spis_XP_022800629_1, Spis_XP_022801342_1, Spis_XP_022801588_1, Spis_XP_022802835_1, Spis_XP_022802839_1, Spis_XP_022802910_1, Spis_XP_022803016_1, Spis_XP_022803357_1, Spis_XP_022804047_1, Spis_XP_022804328_1, Spis_XP_022805991_1, Spis_XP_022806486_1, Spis_XP_022806531_1, Spis_XP_022807308_1, Spis_XP_022807326_1, Spis_XP_022807438_1, Spis_XP_022807804_1, Spis_XP_022808049_1, Spis_XP_022808223_1, Spis_XP_022808237_1, Spis_XP_022808367_1, Spis_XP_022808602_1, Spis_XP_022808640_1, Spis_XP_022808778_1, Spis_XP_022808781_1, Spis_XP_022808881_1, Spis_XP_022808896_1, Spis_XP_022809255_1, Spis_XP_022809256_1, Spis_XP_022809304_1, Spis_XP_022809368_1, Spis_XP_022809407_1, Spis_XP_022809408_1, Spis_XP_022809566_1, Spis_XP_022809651_1, Spis_XP_022809728_1, Spis_XP_022809966_1, Spis_XP_022810044_1, Spis_XP_022810382_1, Spis_XP_022810383_1, Spis_XP_022810486_1, Spis_XP_022810490_1"
OG0000003	transcript_HQ_P2_transcript1460_f2p0_4473_1	"Spis10850_1, Spis11403_1, Spis11947_1, Spis13568_1, Spis1611_1, Spis16461_1, Spis16522_1, Spis16772_1, Spis16973_1, Spis17407_1, Spis18128_1, Spis19094_1, Spis19207_1, Spis19862_1, Spis20535_1, Spis20796_1, Spis21014_1, Spis2122_1, Spis21240_1, Spis22160_1, Spis22271_1, Spis22832_1, Spis23516_1, Spis23661_1, Spis23710_1, Spis23821_1, Spis23912_1, Spis24186_1, Spis24195_1, Spis24208_1, Spis24688_1, Spis25600_1, Spis263_1, Spis307_1, Spis3244_1, Spis3451_1, Spis3742_1, Spis42_1, Spis4580_1, Spis4871_1, Spis5190_1, Spis628_1, Spis6559_1, Spis6681_1, Spis6684_1, Spis6724_1, Spis683_1, Spis8570_1, Spis9042_1, Spis9235_1, Spis9560_1, Spis_XP_022777852_1, Spis_XP_022777960_1, Spis_XP_022778341_1, Spis_XP_022778993_1, Spis_XP_022779436_1, Spis_XP_022779767_1, Spis_XP_022779801_1, Spis_XP_022781137_1, Spis_XP_022782754_1, Spis_XP_022783198_1, Spis_XP_022785752_1, Spis_XP_022785990_1, Spis_XP_022786150_1, Spis_XP_022787108_1, Spis_XP_022788162_1, Spis_XP_022789646_1, Spis_XP_022790920_1, Spis_XP_022792746_1, Spis_XP_022793604_1, Spis_XP_022793726_1, Spis_XP_022793846_1, Spis_XP_022794398_1, Spis_XP_022794661_1, Spis_XP_022796121_1, Spis_XP_022796579_1, Spis_XP_022797031_1, Spis_XP_022797150_1, Spis_XP_022797523_1, Spis_XP_022798804_1, Spis_XP_022798951_1, Spis_XP_022799704_1, Spis_XP_022799813_1, Spis_XP_022800314_1, Spis_XP_022800929_1, Spis_XP_022800936_1, Spis_XP_022802273_1, Spis_XP_022802534_1, Spis_XP_022803076_1, Spis_XP_022803077_1, Spis_XP_022805422_1, Spis_XP_022806191_1, Spis_XP_022806614_1, Spis_XP_022807381_1, Spis_XP_022808390_1, Spis_XP_022808393_1, Spis_XP_022808716_1, Spis_XP_022808790_1, Spis_XP_022809152_1, Spis_XP_022809203_1, Spis_XP_022809437_1, Spis_XP_022810492_1, Spis_XP_022810493_1, Spis_XP_022810656_1"
OG0000004	transcript_HQ_P2_transcript1921_f7p0_4136_1	"Spis11408_1, Spis11697_1, Spis15454_1, Spis17486_1, Spis17720_1, Spis21554_1, Spis22219_1, Spis22807_1, Spis24631_1, Spis_XP_022777497_1, Spis_XP_022777498_1, Spis_XP_022777541_1, Spis_XP_022781085_1, Spis_XP_022781088_1, Spis_XP_022781129_1, Spis_XP_022781133_1, Spis_XP_022781158_1, Spis_XP_022785107_1, Spis_XP_022785539_1, Spis_XP_022792749_1, Spis_XP_022792750_1, Spis_XP_022793079_1, Spis_XP_022793270_1, Spis_XP_022793271_1, Spis_XP_022793273_1, Spis_XP_022793282_1, Spis_XP_022793287_1, Spis_XP_022793667_1, Spis_XP_022793671_1, Spis_XP_022793674_1, Spis_XP_022793683_1, Spis_XP_022794330_1, Spis_XP_022795196_1, Spis_XP_022795733_1, Spis_XP_022796384_1, Spis_XP_022796395_1, Spis_XP_022796405_1, Spis_XP_022796406_1, Spis_XP_022796421_1, Spis_XP_022797339_1, Spis_XP_022797340_1, Spis_XP_022797342_1, Spis_XP_022798301_1, Spis_XP_022800040_1, Spis_XP_022800045_1, Spis_XP_022800750_1, Spis_XP_022801670_1, Spis_XP_022802381_1, Spis_XP_022802382_1, Spis_XP_022802383_1, Spis_XP_022802386_1, Spis_XP_022802389_1, Spis_XP_022802393_1, Spis_XP_022802395_1, Spis_XP_022802710_1, Spis_XP_022802714_1, Spis_XP_022803123_1, Spis_XP_022803124_1, Spis_XP_022803128_1, Spis_XP_022803635_1, Spis_XP_022803659_1, Spis_XP_022803821_1, Spis_XP_022803822_1, Spis_XP_022803836_1, Spis_XP_022804170_1, Spis_XP_022804175_1, Spis_XP_022804176_1, Spis_XP_022804177_1, Spis_XP_022804729_1, Spis_XP_022805120_1, Spis_XP_022805177_1, Spis_XP_022805178_1, Spis_XP_022805604_1, Spis_XP_022805605_1, Spis_XP_022805954_1, Spis_XP_022807724_1, Spis_XP_022807771_1, Spis_XP_022807772_1, Spis_XP_022807774_1, Spis_XP_022807775_1, Spis_XP_022807950_1, Spis_XP_022808095_1, Spis_XP_022808202_1, Spis_XP_022808789_1, Spis_XP_022809167_1, Spis_XP_022809179_1, Spis_XP_022809433_1, Spis_XP_022809434_1, Spis_XP_022809523_1, Spis_XP_022809733_1, Spis_XP_022810254_1, Spis_XP_022810255_1, Spis_XP_022810275_1, Spis_XP_022810324_1, Spis_XP_022810339_1, Spis_XP_022810419_1, Spis_XP_022810529_1"
OG0000005	transcript_HQ_P2_transcript12300_f6p0_2179_1	"Spis12259_1, Spis14562_1, Spis15079_1, Spis15586_1, Spis15808_1, Spis17195_1, Spis18298_1, Spis20281_1, Spis20400_1, Spis20580_1, Spis21288_1, Spis21792_1, Spis21831_1, Spis22272_1, Spis3038_1, Spis321_1, Spis5580_1, Spis7504_1, Spis8347_1, Spis_XP_022777686_1, Spis_XP_022777967_1, Spis_XP_022778516_1, Spis_XP_022780185_1, Spis_XP_022781620_1, Spis_XP_022782206_1, Spis_XP_022782216_1, Spis_XP_022782695_1, Spis_XP_022783276_1, Spis_XP_022783278_1, Spis_XP_022783673_1, Spis_XP_022783684_1, Spis_XP_022783871_1, Spis_XP_022784026_1, Spis_XP_022784269_1, Spis_XP_022784914_1, Spis_XP_022785543_1, Spis_XP_022786864_1, Spis_XP_022786975_1, Spis_XP_022787710_1, Spis_XP_022788757_1, Spis_XP_022789239_1, Spis_XP_022789703_1, Spis_XP_022790235_1, Spis_XP_022790387_1, Spis_XP_022790680_1, Spis_XP_022791336_1, Spis_XP_022791350_1, Spis_XP_022791951_1, Spis_XP_022792007_1, Spis_XP_022792607_1, Spis_XP_022793346_1, Spis_XP_022794086_1, Spis_XP_022795746_1, Spis_XP_022796281_1, Spis_XP_022797384_1, Spis_XP_022797549_1, Spis_XP_022798322_1, Spis_XP_022800734_1, Spis_XP_022800926_1, Spis_XP_022801520_1, Spis_XP_022802259_1, Spis_XP_022802642_1, Spis_XP_022803023_1, Spis_XP_022803257_1, Spis_XP_022803564_1, Spis_XP_022803739_1, Spis_XP_022804685_1, Spis_XP_022804789_1, Spis_XP_022804858_1, Spis_XP_022804901_1, Spis_XP_022805635_1, Spis_XP_022805636_1, Spis_XP_022805897_1, Spis_XP_022806310_1, Spis_XP_022806397_1, Spis_XP_022806610_1, Spis_XP_022806915_1, Spis_XP_022807037_1, Spis_XP_022807099_1, Spis_XP_022807457_1, Spis_XP_022807818_1, Spis_XP_022808002_1, Spis_XP_022808164_1, Spis_XP_022808340_1, Spis_XP_022808560_1, Spis_XP_022808753_1, Spis_XP_022808849_1, Spis_XP_022809248_1, Spis_XP_022810033_1, Spis_XP_022810099_1, Spis_XP_022810141_1, Spis_XP_022810147_1, Spis_XP_022810200_1"
以\t 为连接，每个物种中以", "连接，需要找到上文中的基因对应的1:1的FL物种的基因
"""
def read_orthofinder_results(ortho_file):
    """读取OrthoFinder结果文件，返回1:1对应关系的字典"""
    orthologs_dict = {}
    with open(ortho_file, 'r') as f:
        # 跳过表头
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            # 确保行至少有3个字段
            if len(fields) < 3:
                continue
                
            og_id = fields[0]
            # 处理空字段的情况
            fl_genes = fields[1].split(', ') if fields[1].strip() else []
            sp_genes = fields[2].replace('"', '').split(', ') if fields[2].strip() else []
            
            # 只保留1:1对应关系
            if len(fl_genes) == 1 and len(sp_genes) == 1:
                fl_gene = fl_genes[0]
                sp_gene = sp_genes[0]
                orthologs_dict[sp_gene] = fl_gene
    return orthologs_dict
"""
读取cluster0_marker.tsv文件,读取,只保留p_val < 0.05 且avg_log2FC > 0 的全部基因
"""
def read_cluster_markers(cluster_file):
    """读取cluster marker文件，返回显著性基因列表"""
    significant_genes = []
    with open(cluster_file, 'r') as f:
        # 跳过表头
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            p_val = float(fields[1])
            avg_log2FC = float(fields[2])
            if p_val < 0.05 and avg_log2FC > 0:
                                # 格式化基因名：添加_1后缀，替换/为_，替换-为_
                gene_name = fields[0]
                gene_name = gene_name.replace('/', '_').replace('-', '_')
                if not gene_name.endswith('_1'):
                    gene_name = f"{gene_name}_1"
                significant_genes.append(gene_name)
    return significant_genes


def find_orthologs(mapping_proteins_dict, orthologs_dict, significant_genes):
    """找到对应的1:1直系同源基因"""
    results = []
    for protein, source_files in mapping_proteins_dict.items():
        # 为protein添加Spis_前缀
        spis_protein = f"Spis_{protein}"
        if spis_protein in orthologs_dict:
            fl_gene = orthologs_dict[spis_protein]
            if fl_gene in significant_genes:
                results.append((spis_protein, fl_gene, source_files))
    return results


def main():
    mapping_files = [

        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster0_markers_vs_unknown-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster1_markers_vs_unknown-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster2_markers_vs_unknown-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster3_markers_vs_unknown-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster6_markers_vs_unknown-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_germline-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_gland-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_immune-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_neuron-markers.tsv",
        r"D:\nextcloud\pd论文\result\comparison-result\cluster7_markers_vs_unknown-markers.tsv",
    ]
    ortho_file = r"D:\nextcloud\pd论文\result\Orthofinder-result\FL-SP\FL-SP-cell_Orthogroups.tsv"
    cluster_file = r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv"
    
    # 添加日志输出
    mapping_proteins_dict = read_mapping_file(mapping_files)
    print(f"找到 {len(mapping_proteins_dict)} 个 overlapping proteins")
    
    orthologs_dict = read_orthofinder_results(ortho_file)
    print(f"找到 {len(orthologs_dict)} 个1:1同源基因对")
    
    significant_genes = read_cluster_markers(cluster_file)
    print(f"找到 {len(significant_genes)} 个显著性基因")
    
    results = find_orthologs(mapping_proteins_dict, orthologs_dict, significant_genes)
    print(f"最终找到 {len(results)} 个符合条件的基因对")
    
    # 输出结果
    print("\nProtein ID\tFL Gene\tSource Files")
    print("-" * 60)
    for sp_gene, fl_gene, source_files in results:
        source_str = ", ".join(source_files)  # 将来源文件列表转换为字符串
        print(f"{sp_gene}\t{fl_gene}\t{source_str}")

if __name__ == "__main__":
    main()