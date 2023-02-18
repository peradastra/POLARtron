from opentrons import types
from datetime import datetime
import subprocess
import os
import os.path

metadata = {
    'protocolName': 'POLARtron: Nucleic Acid Extraction & Split Pool RT-PCR Modules',
    'author': 'Per Adastra <adastra.aspera.per@gmail.com>',
    'apiLevel': '2.10'
}

def run(ptx, experiment_name=""):
    run_log_directory = "/var/lib/jupyter/notebooks/run_logs"

    # <editor-fold desc="Create sample dictionary">
    # Dynamically create sample list for run based on the number of samples.
    polar = dict()
    samples = []
    for column in range(1, 5):
        polar['Sample #' + str(column)] = {}
        samples.append('Sample #' + str(column))

    # </editor-fold>

    # <editor-fold desc="Define tips">

    # p200 dynamic tip box alloc
    p200TipRack = [ptx.load_labware('opentrons_96_tiprack_300ul', s) for s in ['9', '6', '5', '8']]

    # p200 tip enumeration
    p200Tips = [tr['A' + str(i)] for tr in p200TipRack for i in reversed(range(1, 13))]

    # p20 tip box alloc
    p20TipRack = ptx.load_labware('opentrons_96_tiprack_20ul', '2')

    # p20 tip enumeration
    p20Tips = [p20TipRack['A' + str(i)] for i in range(1, 13)]

    # </editor-fold>

    # <editor-fold desc="Define hardware">

    # Load instruments
    p300 = ptx.load_instrument('p300_multi_gen2', 'left')
    p20 = ptx.load_instrument('p20_multi_gen2', 'right')

    # Trash bin.
    trash = ptx.load_labware('agilent_1_reservoir_290ml', 11)['A1']

    # Load modules
    magneticModule = ptx.load_module('magnetic module gen2', '4')
    temperatureModule = ptx.load_module('temperature module gen2', '1')
    thermocyclerModule = ptx.load_module('thermocycler')

    # </editor-fold>

    # <editor-fold desc="Define plates, reservoirs and trash">

    # Plates & reservoirs
    reagentResevoir = ptx.load_labware('nest_12_reservoir_15ml', '3')
    thermocyclerPlate = thermocyclerModule.load_labware('biorad_96_wellplate_200ul_pcr')
    coldReagentsPlate = temperatureModule.load_labware('biorad_96_wellplate_200ul_pcr')
    magneticModulePlate = magneticModule.load_labware('eppendorf_96_well_lobind_plate_500ul')

    # </editor-fold>

    # <editor-fold desc="Tip assignments">

    # List of p200 tips needed for each s
    p200TipNeeds = [
        'mixTip',  # Tip 1 (Box 1)
        'viralBufferTip1',  # Tip 2 (Box 1)
        'viralBufferTip2',  # Tip 3 (Box 1)
        'magbeadBufferTip1',  # Tip 4 (Box 2)
        'magbeadBufferTip2',  # Tip 5 (Box 2)
        'ethanolTip1',  # Tip 6 (Box 2)
        'ethanolTip2',  # Tip 7 (Box 3)
        'ethanolTip3',  # Tip 8 (Box 3)
        'elutionTip',  # Tip 9 (Box 3)
        'stopTip',  # Tip 10 (Box 4)
        'bltBeadsWashTip',  # Tip 11 (Box 4)
        'pcrTip'  # Tip 12 (Box 4)
    ]

    # List of p20 tips needed for each s
    p20TipNeeds = [
        'rtPcrPool1Tip',  # Tip 1
        'rtPcrPool2Tip',  # Tip 2
        'bltBeadTip'  # Tip 3
    ]

    # Dynamic tip assignment based on needs listed above
    tipCounter = 0
    for tip in p20TipNeeds:
        for sample in polar:
            polar[sample][tip] = p20Tips[tipCounter]
            tipCounter += 1

    tipCounter = 0
    for tip in p200TipNeeds:
        for sample in polar:
            polar[sample][tip] = p200Tips[tipCounter]
            tipCounter += 1

    tipForMixingAccukitProtinaseK = polar[samples[0]]['viralBufferTip1']
    tipForMixingRtPcr = polar[samples[0]]['ethanolTip1']
    tipForMixingViralBuffer = polar[samples[0]]['viralBufferTip1']
    tipForMixingHackflex = polar[samples[0]]['stopTip']
    tipForAddStopBuffer = polar[samples[0]]['stopTip']

    # </editor-fold>

    # <editor-fold desc="Reagents">

    # room temperature reagents
    viralBufferBeads = reagentResevoir['A1']
    viralBuffer = reagentResevoir['A2']
    magbeadBuffer1 = reagentResevoir['A3']
    magbeadBuffer2 = reagentResevoir['A4']
    ethanolWells = reagentResevoir['A5'], reagentResevoir['A6']
    elutionBuffer = reagentResevoir['A7']
    mineralOil = reagentResevoir['A8']
    stopBuffer = reagentResevoir['A9']
    bltBeadWashBuffer = reagentResevoir['A10']
    wash_well = reagentResevoir['A11']

    index_list = []
    for num in range(9, 13):
        index_list.append('A' + str(num))

    index_counter = 0
    for sample in samples:
        polar[sample]['pcrWithIndex'] = coldReagentsPlate[index_list[index_counter]]
        index_counter += 1

    # </editor-fold>

    # <editor-fold desc="Protocol variables">

    # volume
    rtpcrVolume = 20
    hackflexVolume = 7.5

    # time
    hacklexBeadBindTiime = 5

    # offsets/adjusts/def
    wellFillRate = 0.032
    p20DefaultRate = 7.56
    p300DefaultRate = 92.86
    lobindEngageHeight = 7.4

    # </editor-fold>

    # <editor-fold desc="Well assignments">

    # Magnetic module plate wells
    # Create empty list for extraction, library preparation and tip washing
    magModList = []
    for num in range(1, 13):
        magModList.append('A' + str(num))

    wellCounter = 0
    lobindWellNeeds = ['extractionWell',
                       'bltBeadxWashWell',
                       'sampleTipWashWell']

    for sample in samples:
        for well in lobindWellNeeds:
            polar[sample][well] = magneticModulePlate[magModList[wellCounter]]
            wellCounter += 1

    # Thermocycler plate wells

    # Thermocycler plate wells
    # Create empty list for RT PCR
    rtPcrWells = []
    for num in (3, 4, 5, 6, 7, 8, 9, 10):
        rtPcrWells.append('A' + str(num))

    wellCounter = 0
    thermocyclerWellNeeds = ['rtPcrPool1',
                             'rtPcrPool2']

    for sample in samples:
        for well in thermocyclerWellNeeds:
            polar[sample][well] = thermocyclerPlate[rtPcrWells[wellCounter]]
            wellCounter += 1

    # Create empty list for final PCR
    pcrWells = []
    for num in [1, 2, 11, 12]:
        pcrWells.append('A' + str(num))

    wellCounter = 0
    thermocyclerWellNeeds = ['indexPcrWell']

    for sample in samples:
        for well in thermocyclerWellNeeds:
            polar[sample][well] = thermocyclerPlate[pcrWells[wellCounter]]
            wellCounter += 1

    # </editor-fold>

    # <editor-fold desc="Protocol functions">

    def wash_tip(instrament, wash_well, volume):
        if instrament.current_volume > 0:
            instrament.blow_out(wash_well.bottom(5))
        set_speeds(instrament)
        instrament.mix(5, volume, wash_well.bottom())
        ptx.max_speeds['Z'] = ptx.max_speeds['A'] = 10
        instrament.blow_out(well.top(-10))
        ptx.max_speeds['Z'] = ptx.max_speeds['A'] = None

    def make_directory(directory_path):
        if not os.path.exists(directory_path):
            try:
                os.makedirs(directory_path)
            except Exception:
                pass

    def return_protocol_object(protocol_object):
        protocolObjectString = str(protocol_object).split()[3][:-1]
        return protocolObjectString

    def play_alert_sound(sound="alert"):
        if not ptx.is_simulating():
            soundsPath = 'sounds/'
            # packagePath = return_path_to_package_data()
            playlist = {
                "alert": "alert_sound_1.mp3",
                "stop" : "alert_sound_2.mp3"
            }
            # soundsPath = packagePath + soundsPath + playlist[sound]
            soundsPath = '/var/lib/jupyter/notebooks/' + soundsPath + playlist[sound]
            subprocess.run(['mpg123', soundsPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def update_log(update="", experiment=experiment_name):
        protocolString = return_protocol_object(ptx)
        make_directory(run_log_directory)
        ptx.comment(update)
        if not ptx.is_simulating():
            when = datetime.now().strftime("%y_%m_%d_%H_%M_%S")
            if not experiment == "":
                experiment = experiment + "_"
            logName = experiment + "run_log_" + "_" + protocolString + ".txt"
            runLog = open("/var/lib/jupyter/notebooks/run_logs/" + logName, "a")
            runLog.write(when + ' > ' + update + '\n')
            runLog.close()

    def trash_tip():
        if p300.has_tip:
            pipette = p300
        if p20.has_tip:
            pipette = 20
        set_speeds(pipette)
        if pipette.current_volume > 0:
            pipette.dispense(pipette.current_volume, trash.bottom(5))
            pipette.blow_out()
            slow_exit(pipette, trash)

        pipette.return_tip()

    def set_speeds(pipette, aspirate=0, dispense=0):
        if dispense == 0:
            if pipette == p300:
                dispense = p300DefaultRate
            if pipette == p20:
                dispense = p20DefaultRate
        if aspirate == 0:
            if pipette == p300:
                aspirate = p300DefaultRate
            if pipette == p20:
                aspirate = p20DefaultRate
        pipette.flow_rate.aspirate = aspirate
        pipette.flow_rate.dispense = dispense

    def aspirate_fluid(pipette, volume, location, height=1):
        pipette.aspirate(volume, location.bottom(height))
        ptx.delay(seconds=1)

    def slow_exit(pipette, location, height=0):
        ptx.max_speeds['Z'] = ptx.max_speeds['A'] = 10
        pipette.move_to(location.top(height))
        ptx.max_speeds['Z'] = ptx.max_speeds['A'] = None

    def bead_side(well):
        column = int("".join(filter(str.isdigit, str(well).split()[0])))
        if (column % 2) == 0:
            return -1
        else:
            return 1

    def remove_supernatant(pipette, volume, location):
        side = bead_side(location)
        pipette.flow_rate.aspirate = 10
        pipette.aspirate((volume * 0.7), location.bottom().move(
            types.Point(z=((0.3 * volume) * wellFillRate), x=(-1.5 * side))))
        ptx.delay(seconds=5)
        pipette.flow_rate.aspirate = 10
        pipette.aspirate((volume * 0.30), location.bottom().move(types.Point(z=1, x=(-1.5 * side))))
        set_speeds(p300)

    def resuspend_beads(pipette, reps, volume, location, aspirate=400, dispense=400):
        side = bead_side(location)
        pipette.flow_rate.aspirate = aspirate
        pipette.flow_rate.dispense = dispense
        for i in range(reps):
            pipette.aspirate((volume * 0.8), location.bottom(1))
            pipette.dispense((volume * 0.8), location.bottom(3).move(types.Point(x=1.5 * side)))

    def side_dispense(pipette, location, volume=0, dispense=50, aspirate=50, blowOut=True):
        set_speeds(pipette, aspirate, dispense)
        if volume == 0:
            volume = pipette.current_volume
        pipette.move_to(location.top())
        pipette.default_speed = 25
        pipette.dispense(volume, location.top().move(types.Point(y=4.5, z=-5)))
        ptx.delay(seconds=1)
        if blowOut:
            pipette.blow_out(location.top().move(types.Point(y=4.5, z=-5)))
            ptx.delay(seconds=1)
        pipette.move_to(location.top())
        pipette.default_speed = None
        set_speeds(pipette)

    def well_wash(pipette, location):
        volume = pipette.current_volume / 4
        pipette.flow_rate.aspirate = 50
        pipette.flow_rate.dispense = 50
        pipette.move_to(location.top(-3))
        pipette.default_speed = 20
        for side in (1, -1):
            pipette.move_to(location.top().move(types.Point(x=(side * 3), z=-3)))
            pipette.dispense(volume, location.top().move(types.Point(x=(side * 3), z=-3)))
            pipette.move_to(location.top().move(types.Point(y=(side * 3), z=-3)))
            pipette.dispense(volume, location.top().move(types.Point(y=(side * 3), z=-3)))
        pipette.blow_out()
        ptx.delay(seconds=1)
        pipette.move_to(location.top().move(types.Point(z=0)))
        pipette.default_speed = None

    def engage_magnet_module(minutes=0):
        magneticModule.engage(lobindEngageHeight)
        ptx.delay(minutes=minutes)

    def collect_dispense_touch(pipette, volume, location, aspirate=0, dispense=0,
                               blow_out=False, touch_tip=True):
        set_speeds(pipette, aspirate, dispense)
        pipette.aspirate(volume, location.bottom())
        ptx.delay(seconds=1)
        pipette.dispense(pipette.current_volume, location.bottom(wellFillRate * volume))
        ptx.delay(seconds=1)
        if blow_out:
            slow_exit(pipette, location, height=-10)
            pipette.flow_rate.blow_out = 10
            pipette.blow_out(location.top().move(types.Point(z=-10)))
            ptx.delay(seconds=5)
            pipette.default_speed = None
        if touch_tip:
            well_touch_tip(pipette, location)
        slow_exit(pipette, location)
        set_speeds(pipette)

    def well_touch_tip(pipette, location, height=-5):
        slow_exit(pipette, location, height=height)
        pipette.default_speed = 20
        for side in (1, -1):
            pipette.move_to(location.top().move(types.Point(z=height, x=(side * 4.5))))
            pipette.move_to(location.top().move(types.Point(z=height, y=(side * 4.5))))
        pipette.move_to(location.top(height))
        pipette.default_speed = None

    def wash_beads(pipette, volume, buffer, well, tip, reps=15, time=5, is_detergent=False, resuspend=True):

        if not resuspend:
            reps = 0
            time = 0

        if is_detergent or resuspend:
            aspirate = 100
            dispense = 5
        else:
            aspirate = 0
            dispense = 5

        if resuspend:
            magneticModule.disengage()

        for sample in samples:
            pipette.pick_up_tip(polar[sample][tip])
            aspirate_fluid(pipette, volume, buffer)
            if buffer in ethanolWells:
                well_wash(pipette, polar[sample][well])
            elif not resuspend:
                set_speeds(pipette, 100, 10)
                pipette.dispense(pipette.current_volume, polar[sample][well].bottom())
                slow_exit(pipette, polar[sample][well])
            else:
                side_dispense(pipette, polar[sample][well], dispense=dispense)
                slow_exit(pipette, polar[sample][well])
            pipette.return_tip()

        if resuspend:
            for sample in samples:
                pipette.pick_up_tip(polar[sample][tip])
                resuspend_beads(pipette, reps, volume, polar[sample][well])
                collect_dispense_touch(pipette, volume, polar[sample][well], aspirate=aspirate, dispense=dispense,
                                       blow_out=True)
                pipette.return_tip()

        engage_magnet_module(time)

        for sample in samples:
            pipette.pick_up_tip(polar[sample][tip])
            remove_supernatant(pipette, volume, polar[sample][well])
            slow_exit(p300, polar[sample][well])
            trash_tip()

    def pause_protocol(comment="", sound='default', play_sound=True, required_stop=False):
        if required_stop:
            if play_sound:
                if sound == 'default':
                    play_alert_sound()
                else:
                    play_alert_sound(sound)
            ptx.pause(comment)

    def liquid_level(well_volume):
        height = well_volume * wellFillRate
        return height

    # </editor-fold>

    """
    Protocol starts below.
    """

    # <editor-fold desc="Set up OT2 and modules for run.">
    update_log("ʕ·ᴥ·ʔ : OT-2 module set up started.")

    ptx.set_rail_lights(True)
    thermocyclerModule.open_lid()
    engage_magnet_module()
    update_log("ʕ·ᴥ·ʔ : Cooling thermocycler plate to 4°C.")
    thermocyclerModule.set_block_temperature(4)
    update_log("ʕ·ᴥ·ʔ : Cooling temperature module to 4°C.")
    temperatureModule.set_temperature(4)

    update_log("ʕ·ᴥ·ʔ : OT-2 module set up complete.")
    # </editor-fold> #

    # <editor-fold desc="Place samples onto OT-2.">
    update_log("ʕ·ᴥ·ʔ : Awaiting samples to be loaded.")

    pause_protocol("Place sample plate onto magnetic module and press resume to begin.", required_stop=True)
    magneticModule.disengage()

    update_log("ʕ·ᴥ·ʔ : Samples loaded, protocol started.")
    # </editor-fold> #

    # <editor-fold desc="Add extraction control and Protinase K.">
    update_log("ʕ·ᴥ·ʔ : Adding extraction control and Protinase K.")

    update_log("ʕ·ᴥ·ʔ : Mixing Protinase K & Accukit master mix.")
    p300.pick_up_tip(tipForMixingAccukitProtinaseK)
    set_speeds(p300)
    p300.mix(30, 100, coldReagentsPlate['A1'].bottom(1.5))
    slow_exit(p300, coldReagentsPlate['A1'], -2.5)
    p300.flow_rate.blow_out = 10
    p300.blow_out()
    well_touch_tip(p300, coldReagentsPlate['A1'], -2.5)
    p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Adding Protinase K & Accukit master mix to each sample.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['mixTip'])
        set_speeds(p300)
        aspirate_fluid(p300, 25, coldReagentsPlate['A1'])
        slow_exit(p300, coldReagentsPlate['A1'])
        p300.dispense(p300.current_volume, polar[sample]['extractionWell'].bottom())
        set_speeds(p300, 400, 400)
        p300.mix(30, 90)
        collect_dispense_touch(p300, 90, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Incubating sample with Protinase K.")
    ptx.delay(minutes=10)

    update_log("ʕ·ᴥ·ʔ : Extraction control added and Protinase K treatment complete")
    # </editor-fold> #

    # <editor-fold desc="Plate RT-PCR reactions.">
    update_log("ʕ·ᴥ·ʔ : Plating RT-PCR reactions.")

    update_log("ʕ·ᴥ·ʔ : Mixing RT-PCR master mixes.")
    p300.pick_up_tip(tipForMixingRtPcr)
    set_speeds(p300)
    p300.mix(30, 80, coldReagentsPlate['A4'].bottom(1.5))
    slow_exit(p300, coldReagentsPlate['A4'], -2.5)
    p300.flow_rate.blow_out = 10
    p300.blow_out()
    well_touch_tip(p300, coldReagentsPlate['A4'], -2.5)
    p300.return_tip()

    for pool in ('rtPcrPool1', 'rtPcrPool2'):
        if pool == 'rtPcrPool1':
            msg = "ʕ·ᴥ·ʔ : Plating Pool 1 RT-PCR master mix."
            rtPcrPoolMix = coldReagentsPlate['A4']
            p20.pick_up_tip(p20TipRack['E4'])
        if pool == 'rtPcrPool2':
            msg = "ʕ·ᴥ·ʔ : Plating Pool 2 RT-PCR master mix."
            rtPcrPoolMix = coldReagentsPlate['E4']
            p20.pick_up_tip(p20TipRack['E5'])

        update_log(msg)
        for sample in samples:
            aspirate_fluid(p20, 12.5, rtPcrPoolMix, height=0.5)
            slow_exit(p20, rtPcrPoolMix)
            p20.dispense(12.5, polar[sample][pool].bottom(1).move(types.Point(y=-36)))
            ptx.max_speeds['Z'] = ptx.max_speeds['A'] = 10
            p20.move_to(polar[sample][pool].top().move(types.Point(y=-36)))
            ptx.max_speeds['Z'] = ptx.max_speeds['A'] = None
            aspirate_fluid(p20, 12.5, rtPcrPoolMix, height=0.5)
            slow_exit(p20, rtPcrPoolMix)
            p20.move_to(polar[sample][pool].top())
            p20.dispense(12.5, polar[sample][pool].bottom())
            slow_exit(p20, polar[sample][pool])

        p20.return_tip()

    update_log("ʕ·ᴥ·ʔ : Closing thermocycler lid and deactivating temperature module.")
    temperatureModule.deactivate()
    thermocyclerModule.close_lid()
    ptx.home()

    update_log("ʕ·ᴥ·ʔ : RT-PCR reaction plated.")
    # </editor-fold>

    # <editor-fold desc="Binding DNA/RNA to MagBeads.">
    update_log("ʕ·ᴥ·ʔ : Bind DNA/RNA to MagBeads.")

    update_log("ʕ·ᴥ·ʔ : Resuspending MagBeads in Viral DNA/RNA Buffer.")
    p300.pick_up_tip(tipForMixingViralBuffer)
    set_speeds(p300, 400, 400)
    for _ in range(60):
        p300.aspirate(180, viralBufferBeads.bottom())
        p300.dispense(p300.current_volume, viralBufferBeads.bottom(5))
    slow_exit(p300, viralBufferBeads)  # TODO Add blow out and touch wall of well.
    p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Adding Viral DNA/RNA Buffer with MagBeads to all samples..")
    for sample in samples:
        set_speeds(p300, 100, 5)
        p300.pick_up_tip(polar[sample]['viralBufferTip1'])
        aspirate_fluid(p300, 125, viralBufferBeads)
        slow_exit(p300, viralBufferBeads)
        p300.dispense(p300.current_volume, polar[sample]['extractionWell'].bottom(liquid_level(250)))
        slow_exit(p300, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Adding Viral DNA/RNA Buffer to all samples..")
    for sample in samples:
        set_speeds(p300, 100, 5)
        p300.pick_up_tip(polar[sample]['viralBufferTip2'])
        aspirate_fluid(p300, 125, viralBuffer)
        slow_exit(p300, viralBuffer)
        p300.dispense(p300.current_volume, polar[sample]['extractionWell'].bottom(liquid_level(375)))
        slow_exit(p300, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Mixing Viral DNA/RNA Buffer with MagBeads into all samples.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['viralBufferTip1'])
        set_speeds(p300)
        for _ in range(30):
            p300.aspirate(180, polar[sample]['extractionWell'].bottom())
            p300.dispense(p300.current_volume, polar[sample]['extractionWell'].bottom(5))
        collect_dispense_touch(p300, 180, polar[sample]['extractionWell'], dispense=5, blow_out=True)
        slow_exit(p300, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Allow DNA/RNA to bind to MagBeads.")
    ptx.delay(minutes=10)

    update_log("ʕ·ᴥ·ʔ : Pelleting MagBeads.")
    engage_magnet_module(minutes=12)

    update_log("ʕ·ᴥ·ʔ : Remove Viral DNA/RNA Buffer.")
    for tip in ['viralBufferTip1', 'viralBufferTip2']:
        for sample in samples:
            side = bead_side(polar[sample]['extractionWell'])
            p300.pick_up_tip(polar[sample][tip])
            p300.flow_rate.aspirate = 50
            p300.move_to(polar[sample]['extractionWell'].top())
            p300.aspirate(180, polar[sample]['extractionWell'].bottom().move(
                types.Point(z=1, x=(-2 * side))))
            p300.aspirate(20, polar[sample]['extractionWell'].bottom().move(
                types.Point(z=0.5, x=(-2 * side))))
            ptx.delay(seconds=1)
            slow_exit(p300, polar[sample]['extractionWell'])
            trash_tip()

    update_log("ʕ·ᴥ·ʔ : DNA/RNA bound to MagBeads.")
    # </editor-fold>

    # <editor-fold desc="Wash MagBeads with MagBead Wash Buffers 1 & 2.">
    update_log("ʕ·ᴥ·ʔ : Washing MagBeads with MagBead Wash Buffers 1 & 2.")

    # Wash beads in Magbead Buffer 1 and Magbead Buffer 2.
    magbeadBuffersTips = ('magbeadBufferTip1', 'magbeadBufferTip2')
    magbeadBuffers = (magbeadBuffer1, magbeadBuffer2)
    for tip, buffer in zip(magbeadBuffersTips, magbeadBuffers):
        wash_beads(p300, 150, buffer, 'extractionWell', tip, 20, is_detergent=True)

    update_log("ʕ·ᴥ·ʔ : MagBeads washed with MagBead Wash Buffers 1 & 2.")
    # </editor-fold>

    # <editor-fold desc="Wash MagBeads with ethanol.">
    update_log("ʕ·ᴥ·ʔ : Washing MagBeads with ethanol.")

    # Wash MagBeads with ethanol.
    ethanolTips = ['ethanolTip1', 'ethanolTip2']
    for tip, buffer in zip(ethanolTips, ethanolWells):
        wash_beads(p300, 175, buffer, 'extractionWell', tip, 20)

    update_log("ʕ·ᴥ·ʔ : Removing residual amounts of ethanol left in well.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['ethanolTip3'])
        p300.aspirate(50, polar[sample]['extractionWell'].bottom())
        p300.aspirate(50, polar[sample]['extractionWell'].bottom(-1))
        slow_exit(p300, polar[sample]['extractionWell'], height=-15)
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Allowing MagBeads to dry.")
    magneticModule.disengage()
    ptx.delay(minutes=5)

    update_log("ʕ·ᴥ·ʔ : MagBeads washed with ethanol.")
    # </editor-fold>

    # <editor-fold desc="Elute DNA/RNA from MagBeads.">
    update_log("ʕ·ᴥ·ʔ : Eluting DNA/RNA from MagBeads.")

    update_log("ʕ·ᴥ·ʔ : Opening thermocycler lid.")
    thermocyclerModule.open_lid()

    update_log("ʕ·ᴥ·ʔ : Adding Elution Buffer to all wells.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['elutionTip'])
        aspirate_fluid(p300, 20, elutionBuffer)
        slow_exit(p300, elutionBuffer)
        p300.dispense(p300.current_volume, polar[sample]['extractionWell'].bottom())
        slow_exit(p300, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Mixing MagBeads into Elution Buffer.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['elutionTip'])
        set_speeds(p300, 400, 400)
        p300.mix(20, 16, polar[sample]['extractionWell'].bottom())
        slow_exit(p300, polar[sample]['extractionWell'])
        p300.return_tip()

    update_log("ʕ·ᴥ·ʔ : Pelleting MagBeads.")
    engage_magnet_module(minutes=3)

    update_log("ʕ·ᴥ·ʔ : DNA/RNA eluted from MagBeads.")
    # </editor-fold>

    # <editor-fold desc="Transfer eluent to thermocycler.">
    update_log("ʕ·ᴥ·ʔ : Transfering eluent to thermocycler.")

    # Transfer sample needed for RT-PCR into thermocycler.
    pools = ('rtPcrPool1', 'rtPcrPool2')
    poolTips = ('rtPcrPool1Tip', 'rtPcrPool2Tip')
    for tip, pool in zip(poolTips, pools):
        for sample in samples:
            p20.pick_up_tip(polar[sample][tip])
            side = bead_side(polar[sample]['extractionWell'])
            p20.move_to(polar[sample]['extractionWell'].top(-10))
            p20.flow_rate.aspirate = 5
            p20.aspirate(7.5, polar[sample]['extractionWell'].bottom().move(types.Point(z=0, x=(-3 * side))))
            slow_exit(p20, polar[sample]['extractionWell'])
            p20.dispense(p20.current_volume, polar[sample][pool].bottom())
            set_speeds(p20, 20, 20)
            p20.mix(5, 15, polar[sample][pool].bottom())
            slow_exit(p20, polar[sample][pool])
            p20.return_tip()

    update_log("ʕ·ᴥ·ʔ : Eluent transfer to thermocycler complete.")
    # </editor-fold>

    # <editor-fold desc="Add mineral oil overlay to RT-PCR reactions.">
    update_log("ʕ·ᴥ·ʔ : Adding mineral oil overlay to RT-PCR reactions.")
    for sample in samples:
        p300.pick_up_tip(polar[sample]['bltBeadsWashTip'])
        p300.aspirate(65, mineralOil.bottom())
        slow_exit(p300, mineralOil)
        for pool in pools:
            side_dispense(p300, polar[sample][pool], volume=30, blowOut=False)
            slow_exit(p300, polar[sample][pool])
        p300.return_tip()

    magneticModule.disengage()
    update_log("ʕ·ᴥ·ʔ : Mineral oil overlay added to RT-PCR reactions.")
    # </editor-fold>

    # <editor-fold desc="RT-PCR">
    update_log("ʕ·ᴥ·ʔ : Performing RT-PCR.")

    update_log("ʕ·ᴥ·ʔ : Closing thermocycler lid.")
    thermocyclerModule.close_lid()
    ptx.home()
    thermocyclerModule.set_lid_temperature(105)

    update_log("ʕ·ᴥ·ʔ : Performing uracil DNA glycosylase sample pre-treatment.")
    thermocyclerModule.set_block_temperature(25, block_max_volume=50, hold_time_minutes=3)
    update_log("ʕ·ᴥ·ʔ : Performing reverse transcription.")
    thermocyclerModule.set_block_temperature(55, block_max_volume=50, hold_time_minutes=15)
    update_log("ʕ·ᴥ·ʔ : Performing reverse transcription.")
    thermocyclerModule.set_block_temperature(95, block_max_volume=50, hold_time_minutes=2)
    update_log("ʕ·ᴥ·ʔ : Performing amplicon generation.")
    pcr_profile = [
        {'temperature': 95, 'hold_time_seconds': 15},
        {'temperature': 63, 'hold_time_seconds': 180}]
    thermocyclerModule.execute_profile(steps=pcr_profile, repetitions=30, block_max_volume=50)

    thermocyclerModule.set_block_temperature(4, block_max_volume=50)
    update_log("ʕ·ᴥ·ʔ : RT-PCR complete.")

    # </editor-fold>