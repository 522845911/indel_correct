#!/usr/bin/env python
# coding:utf-8
from __future__ import division
import re, json, sys,Error,BedAnno
sys.dont_write_bytecode = True

class BedAnnoAnno(object):
    def __init__(self, var):
        self.var = var

    @staticmethod
    def TO_JSON(data):
        return json.dumps(data)

    @staticmethod
    def getSOfromR(r):
        if r == '?':
            return 'annotation-fail'
        elif r == 'PROM':
            return '167'
        elif re.match(r"5U", r):
            return '204'
        elif re.match(r"3U", r):
            return '205'
        elif re.match(r"C", r):
            return '316'
        elif re.match(r"R", r):
            return '655'
        elif re.match(r"D", r):
            return '163'
        elif re.match(r"A", r):
            return '164'
        elif re.match(r"I5U", r):
            return '447'
        elif re.match(r"I3U", r):
            return '448'
        elif re.match(r"IC", r) or re.match(r"IR", r):
            return '191'
        else:  # default to intergenic_region
            return '605'

    '''
    =head2 cal_hgvs_pos

    About   : calculate nDot, cDot HGVS position, depend on given offset,
              assign trAlt string and nDot HGVS and cDot HGVS positions.
    Usage   : $annoEnt->cal_hgvs_pos(
                    offset => $offset,
                    tid    => $tid,
                    LR     => $lr,
                    tidDetail => $rh_tidDetail,
              );
    Args    : ofst is total offset to the start(left) of currunt annoblk,
              tid is the transcript id for the detail entry
              tidDetail is the currunt annoblk detail
              LR indicate this offset is left or right pos,
                1 for left and assign sta group,
                0 for right and assign end group.
              "noassign" to indicate whether to assign those information
              to annoEnt, it's used return the useful information only,
              without change the annoEnt.
    Returns : if noassign is used, then return a hash ref, which contains
                { nDot, cDot, exin, r } if successful.
              otherwise, 0 for no assigned status, 1 for successful assigned.

    Notes   : For position mirror on transcript, there are 2 other cases
              than normal case:
              1. annotation fail, which can not be annotated in the region
                 of it, the bad alignment string start with 'E'.
              2. block with length changing mismatches, or long substitution
                 mismatch, which contain the following three cases:

                 a. insertion (I) on refSeq

                        +-------+---+--------+  refgenome
                         \       \ /        /
                          +-------+--------+    refSeq

                 b. deletion (D) on refSeq

                          +-------+--------+    refgenome
                         /       / \        \
                        +-------+---+--------+  refSeq

                 c. delins (S/DI) on refSeq (equal/non-equal length)

                        +-------+---+--------+  refgenome
                        |       |  /        /
                        +-------+-+--------+    refSeq

                 Insertion will have an reversed start/end position on refSeq,
                 due to the 1-based position description system.

                 Any position located in a non-zero length refgenome mismatch
                 block have to extend to total region of mismatched block,
                 and alternate "trAlt" value in annotation entry for this tid.

              This method assign the following tag to $annoEnt
                {
                    trInfo => {
                        $tid => {
                            rnaBegin, rnaEnd,  cdsBegin, cdsEnd,
                            ei_Begin, ei_End,  r_Begin,  r_End,
                            genepartSO, trAlt,
                          },
                        ...
                    }
                }

    =cut
    '''

    def cal_hgvs_pos(annoEnt, tid=None, tidDetail=None, offset=None, LR=None, noassign=None):
        if noassign is not None and int(noassign)==0:
            del noassign
        ofst, tid, rtidDetail, lr = offset, tid, tidDetail, LR
        nDot, cDot = "", ""
        strd = 1 if rtidDetail['strd'] == '+' else 0
        if noassign is None:
            trAlt = annoEnt.trInfo[tid]['trAlt']

        exin = rtidDetail['exin']
        blka = rtidDetail['blka']
        gpSO = rtidDetail['gpSO']

        lofst = ofst
        lofst = int(lofst)
        if lr:  # only one time assignment for one var, offset for left
            rofst = int(rtidDetail['wlen']) - int(ofst) - 1
            if lofst < 0:
                if strd:

                    # 5'promoter is always available
                    nDot = int(rtidDetail['nsta']) + lofst
                    m = re.match(r'(\S+)\-u(\d+)', rtidDetail['csta'])
                    if m:
                        cDot = "{}{}{}".format(m.group(1), '-u', (int(m.group(2)) - lofst))
                else:  # outside transcript region into 3'downstream
                    nDot = "{}{}".format('+', (0 - lofst))
                    m = re.match(r'\*?\d+$', rtidDetail['csta'])
                    if m:
                        cDot = "{}{}{}".format(rtidDetail['csta'], '+d', (0 - lofst))

                    exin = '.'
                    blka = '3D'
                    gpSO = '605'
            elif lofst > int(rtidDetail['wlen']):
                if noassign is None:
                    raise Error.Error(
                        "{}{}{}{}{}".format("Error: exceed the whole length [", lofst, ">", rtidDetail['wlen'], "]"))
                else:
                    return 0
            else:  # 0 <= $lofst <= $$rtidDetail{wlen}

                # 1. check if a mismatch block (only exon have this)
                if rtidDetail['mismatch'] != "" and rtidDetail['mismatch'] != "?":
                    # for each kind of mismatch block,
                    # hit the left edge of var, the refSeq ref will
                    # contain start from the start of the mismatch
                    # block, case I in database already has
                    # a reversed coordinate between left and right.
                    # transcript-alt string shoule be assigned to
                    # this tid
                    nDot = rtidDetail['nsta']
                    cDot = rtidDetail['csta']
                    if noassign is None and lofst > 0:

                        # mismatch info record the transcript-stranded
                        # reference sequence
                        mType, mStart, mEnd, strand_ref = rtidDetail['mismatch'].split(",")
                        if mType == 'E':  # edge cross alignment
                            trAlt = '?'
                        if strd:  # cut left
                            if trAlt != '?':
                                trAlt = "{}{}".format(strand_ref[:lofst], trAlt)
                        else:  # cut right
                            if trAlt != '?':
                                trAlt = "{}{}".format(trAlt, strand_ref[-lofst:])

                # 2. check if hit a normal block's right edge
                elif lofst == int(rtidDetail['wlen']):

                    # give back the chance for left edge parsing
                    if noassign is None:
                        del annoEnt.trInfo[tid]
                    return 0

                # 3. check if annotation-fail
                elif str(gpSO) == 'annotation-fail':  # may not exist from 0.73
                    pass

                # 4. check if a promoter
                elif str(blka) == 'PROM':
                    nDot = int(rtidDetail['nsta']) + int(lofst) if strd else int(rtidDetail['nsta']) - int(lofst)
                    m = re.match(r'(\S+)\-u(\d+)$', rtidDetail['csta'])
                    if m:
                        cDot = m.group(1) + '-u' + str(int(m.group(2)) - lofst) if strd else m.group(1) + '-u' + str(
                        int(m.group(2)) + lofst)

                # 5. check if an exon
                elif re.search(r'EX', exin):
                    nDot = int(rtidDetail['nsta']) + int(lofst) if strd else int(rtidDetail['nsta']) - int(lofst)
                    m = re.match(r'(\*?)(\-?\d+)$', rtidDetail['csta'])
                    if m:
                        cDot = "{}{}".format(m.group(1), (int(m.group(2)) + lofst)) if strd else "{}{}".format(
                            m.group(1), (int(m.group(2)) - lofst))

                # 6. intron
                elif re.search(r'IVS', exin):
                    half_length = int(rtidDetail['wlen']) / 2
                    if str(gpSO) == 'abnormal-intron':

                        # for abnormal intron
                        # the length may be less than 2 bp
                        # and then the nsta nsto and csta csto
                        # will point to the neighbor exons' edge
                        if int(lofst) >= int(half_length):
                            if re.search(r'\d+', rtidDetail['nsto']):
                                nDot = "{}{}{}".format(rtidDetail['nsto'], "-",
                                                       (int(rofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['nsto'], "+", (int(rofst) + 1))
                            if re.search(r'\d+', rtidDetail['csto']):
                                cDot = "{}{}{}".format(rtidDetail['csto'], "-",
                                                       (int(rofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['csto'], "+", (int(rofst) + 1))
                        else:
                            if re.search(r'\d+', rtidDetail['nsta']):
                                nDot = "{}{}{}".format(rtidDetail['nsta'], "+",
                                                       (int(lofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['nsta'], "-", (int(lofst) + 1))
                            if re.search(r'\d+', rtidDetail['csta']):
                                cDot = rtidDetail['csta'] + "+" + (int(lofst) + 1) if strd else rtidDetail[
                                                                                                    'csta'] + "-" + (
                                                                                                int(lofst) + 1)
                    else:
                        if int(lofst) >= int(half_length):  # drop into latter part
                            m = re.match(r'(\d+\+?)(\-?\d+)$', rtidDetail['nsto'])
                            if m:
                                nDot = "{}{}".format(m.group(1), (int(m.group(2)) - rofst)) if strd else "{}{}".format(
                                    m.group(1) , (int(m.group(2)) + rofst))
                            m = re.match(r'([\*\-]?\d+\+?)(\-?\d+)$', rtidDetail['csto'])
                            if m:
                                cDot = "{}{}".format(m.group(1), (int(m.group(2)) - rofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) + rofst))
                        else:
                            m = re.match(r'(\d+\+?)(\-?\d+)$', rtidDetail['nsta'])
                            if m:
                                nDot = "{}{}".format(m.group(1), (int(m.group(2)) + lofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) - lofst))
                            m = re.match(r'([\*\-]?\d+\+?)(\-?\d+)$', rtidDetail['csta'])
                            if m:
                                cDot = "{}{}".format(m.group(1), (int(m.group(2)) + lofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) - lofst))
                else:
                    raise Error.Error("Error: what's this? " + rtidDetail - ['exin'] + "\n")

            if noassign is None:

                # the Begin End assignment is always from 5'->3' except for insertion
                if strd:
                    annoEnt.trInfo[tid]['rnaBegin'] = nDot
                    annoEnt.trInfo[tid]['cdsBegin'] = cDot
                    annoEnt.trInfo[tid]['ei_Begin'] = exin
                    annoEnt.trInfo[tid]['r_Begin'] = blka
                else:
                    annoEnt.trInfo[tid]['rnaEnd'] = nDot
                    annoEnt.trInfo[tid]['cdsEnd'] = cDot
                    annoEnt.trInfo[tid]['ei_End'] = exin
                    annoEnt.trInfo[tid]['r_End'] = blka
                annoEnt.trInfo[tid]['genepartSO'] = gpSO
                annoEnt.trInfo[tid]['trAlt'] = trAlt
        else:  # can be multiple times assignment if large var, offset for right
            real_ofst = int(lofst) - 1  # end_ofst change
            rofst = int(rtidDetail['wlen']) - ofst
            if int(lofst) > int(rtidDetail['wlen']):
                ex_ofst = int(lofst) - int(rtidDetail['wlen'])
                if strd:
                    nDot = "{}{}".format('+', ex_ofst)
                    m = re.match(r'\*?\d+$', rtidDetail['csto'])
                    if m:
                        cDot = "{}{}{}".format(rtidDetail['csto'], "+d", ex_ofst)
                    exin = '.'
                    blka = '3D'
                    gpSO = '605'
                else:  # outside 5'promoter region ?
                    m = re.match(r'\-\d+$', rtidDetail['nsto'])
                    nDot = rtidDetail['nsto'] - ex_ofst if m else -ex_ofst
                    m = re.match(r'(\S+)\-u(\d+)', rtidDetail['csto'])
                    if m:
                        cDot = "{}{}{}".format(m.group(1), '-u', (int(m.group(2)) + ex_ofst))
            elif int(lofst) < 0:  # not proper called
                if noassign is None:
                    raise Error.Error("Error: right preceed the block [" + str(lofst) + " < 0]")
                else:
                    return 0
            else:  # 0 <= $lofst <= $$rtidDetail{wlen}
                # 1. check if a mismatch block (exon)
                if rtidDetail['mismatch'] != "" and rtidDetail['mismatch'] != "?":
                    nDot = rtidDetail['nsto']
                    cDot = rtidDetail['csto']
                    if noassign is None and int(lofst) < int(rtidDetail['wlen']):
                        mType, mStart, mEnd, strand_ref = rtidDetail['mismatch'].split(",")
                        if mType == 'E':  # edge cross alignment
                            trAlt = '?'
                        cut_len = int(rtidDetail['wlen']) - int(lofst)
                        if strd:  # cut right
                            if trAlt != '?':
                                trAlt = "{}{}".format(trAlt, strand_ref[-cut_len:])
                        else:  # cut left
                            if trAlt != '?':
                                trAlt = "{}{}".format(strand_ref[:cut_len], trAlt)

                # 2. check if hit a normal block's left edge
                elif int(lofst) == 0:
                    return 0

                # 3. check if annotation-fail
                elif str(gpSO) == 'annotation-fail':
                    pass  # ( $nDot, $cDot ) = ( '?', '?' )

                # 4. check if a promoter
                elif re.match(r'PROM', blka):
                    nDot = int(rtidDetail['nsta']) + real_ofst if strd else int(rtidDetail['nsta']) - real_ofst
                    m = re.match(r'(\S+)\-u(\d+)$', rtidDetail['csta'])
                    if m:
                        if strd:
                            cDot = "{}{}{}".format(m.group(1), '-u', (int(m.group(2)) - real_ofst))
                        else:
                            cDot = "{}{}{}".format(m.group(1), '-u', (int(m.group(2)) + real_ofst))

                # 5. check if an exon
                elif re.search(r'EX', exin):
                    nDot = int(rtidDetail['nsta']) + real_ofst if strd else int(rtidDetail['nsta']) - real_ofst
                    m = re.match(r'(\*?)(\-?\d+)$', rtidDetail['csta'])
                    if m:
                        cDot = "{}{}".format(m.group(1), (int(m.group(2)) + real_ofst)) if strd else "{}{}".format(
                            m.group(1), (int(m.group(2)) - real_ofst))

                # 6. intron
                elif re.search(r'IVS', exin):
                    half_length = int(rtidDetail['wlen']) / 2
                    if gpSO == 'abnormal-intron':
                        if int(real_ofst) >= int(half_length):
                            m = re.search(r'\d+', rtidDetail['nsto'])
                            if m:
                                nDot = "{}{}{}".format(rtidDetail['nsto'], "-",
                                                       (int(rofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['nsto'], "+", (int(rofst) + 1))
                            m = re.search(r'\d+', rtidDetail['csto'])
                            if m:
                                cDot = "{}{}{}".format(
                                    rtidDetail['csto'] + "-" + (int(rofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['csto'] + "+" + (int(rofst) + 1))
                        else:
                            m = re.search(r'\d+', rtidDetail['nsta'])
                            if m:
                                nDot = "{}{}{}".format(rtidDetail['nsta'], "+",
                                                       (int(rofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['nsta'], "-", (int(rofst) + 1))
                            m = re.search(r'\d+', rtidDetail['csta'])
                            if m:
                                cDot = "{}{}{}".format(rtidDetail['csta'], "+",
                                                       (int(real_ofst) + 1)) if strd else "{}{}{}".format(
                                    rtidDetail['csta'], "-", (int(real_ofst) + 1))
                    else:
                        if int(real_ofst) >= int(half_length):  # drop into latter part
                            m = re.match(r'(\d+\+?)(\-?\d+)$', rtidDetail['nsto'])
                            if m:
                                nDot = "{}{}".format(m.group(1), (int(m.group(2)) - rofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) + rofst))
                            m = re.match(r'([\*\-]?\d+\+?)(\-?\d+)$', rtidDetail['csto'])
                            if m:
                                cDot = "{}{}".format(m.group(1), (int(m.group(2)) - rofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) + rofst))
                        else:
                            m = re.match(r'(\d+\+?)(\-?\d+)$', rtidDetail['nsta'])
                            if m:
                                nDot = "{}{}".format(m.group(1), (int(m.group(2)) + real_ofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) - real_ofst))
                            m = re.match(r'([\*\-]?\d+\+?)(\-?\d+)$', rtidDetail['csta'])
                            if m:
                                cDot = "{}{}".format(m.group(1), (int(m.group(2)) + real_ofst)) if strd else "{}{}".format(
                                    m.group(1), (int(m.group(2)) - real_ofst))
                else:
                    raise Error.Error("Error: what's this? " + rtidDetail["exin"] + "\n")

            if noassign is None:

                # the Begin End assignment is always from 5'->3' except for insertion
                if strd:
                    annoEnt.trInfo[tid]['rnaEnd'] = nDot
                    annoEnt.trInfo[tid]['cdsEnd'] = cDot
                    annoEnt.trInfo[tid]['ei_End'] = exin
                    annoEnt.trInfo[tid]['r_End'] = blka
                else:
                    annoEnt.trInfo[tid]['rnaBegin'] = nDot
                    annoEnt.trInfo[tid]['cdsBegin'] = cDot
                    annoEnt.trInfo[tid]['ei_Begin'] = exin
                    annoEnt.trInfo[tid]['r_Begin'] = blka
                annoEnt.trInfo[tid]['trAlt'] = trAlt
                if rtidDetail['gpSO'] == 'annotation-fail' or (
                        "genepartSO" in annoEnt.trInfo[tid] and annoEnt.trInfo[tid]['genepartSO'] == 'annotation-fail'):
                    annoEnt.trInfo[tid]['genepartSO'] = 'annotation-fail'
                elif "genepartSO" in annoEnt.trInfo[tid] and annoEnt.trInfo[tid]['genepartSO'] != "span" and rtidDetail[
                    'gpSO'] != annoEnt.trInfo[tid]['genepartSO']:
                    annoEnt.trInfo[tid]['genepartSO'] = 'span'
                elif "genepartSO" not in annoEnt.trInfo[tid]:
                    print "Warning: no genepartSO specified in left? line ", sys._getframe().f_lineno
                    annoEnt.trInfo[tid]['genepartSO'] = gpSO
        if noassign is not None:
            return {"nDot": nDot, "cDot": cDot, "exin": exin, "r": blka}
        return 1

    def getTrPosition(annoEnt, rannodb, aeIndex):
        var = annoEnt.var

        # gather the covered entries
        new_aeIndex = None

        tidExblk = dict()
        rpreLeft = dict()

        # debug
        #    print STDERR "DEBUG: getTrPosition for ".Dumper($annoEnt)

        hit_badaln_ins = dict()
        for k in range(aeIndex, len(rannodb)):
            if int(rannodb[k]['sto']) < int(var.pos):  # not reach var
                aeIndex += 1
                continue
            elif int(rannodb[k]['sta']) > int(var.end):  # past var
                break
            else:  # covered by var

                #	    # debug info
                #	    print STDERR "blk: $k [ $$rannodb[$k]{sta}, $$rannodb[$k]{sto} ]\n";

                if int(rannodb[k]['sta']) <= int(var.pos) and int(var.pos) <= int(rannodb[k]['sto']):
                    new_aeIndex = k

                    if new_aeIndex > 0 and int(rannodb[k]['sta']) == int(var.pos):
                        # for edge hit case fetch backward one block.
                        new_aeIndex -= 1
                if "detail" not in rannodb[k]:  # parse anno db
                    rannodb[k] = BedAnno.BedAnno.assign_detail(rannodb[k])

                for tid in sorted(rannodb[k]['detail'].keys()):
                    rtidDetail = rannodb[k]['detail'][tid]
                    strd = None
                    strd = 1 if rtidDetail['strd'] == '+' else 0
                    unify_p, unify_r, unify_a, unify_rl, unify_al = var.getUnifiedVar(rtidDetail['strd'])

                    if int(unify_p) > int(rannodb[k]['sto']) or int(int(unify_p) + int(unify_rl)) < int(rannodb[k]['sta']):
                        # skip non hitted block
                        continue

                    if int(rannodb[k]['sta']) == int(rannodb[k]['sto']) and int(int(unify_p) + int(unify_rl)) == int(rannodb[k]['sto']) and re.search(r'D', rtidDetail['mismatch']):
                        hit_badaln_ins[tid] = 1

                    total_left_ofst = int(rtidDetail['offset']) + int(unify_p) - int(rannodb[k]['sta'])
                    total_right_ofst = int(total_left_ofst) + int(unify_rl)

                    # debug
                    #		print STDERR join(" ", "tid: $tid", $rtidDetail->{strd},
                    #		    $rtidDetail->{nsta}, $rtidDetail->{nsto},
                    #		    $rtidDetail->{mismatch})."\n";
                    #		print STDERR "rel: $total_left_ofst $total_right_ofst\n";

                    if not hasattr(annoEnt,"trInfo") or tid not in annoEnt.trInfo or "trAlt" not in annoEnt.trInfo[tid]:
                        # $tid added into trInfo
                        if not hasattr(annoEnt, "trInfo"):
                            annoEnt.trInfo=dict()
                        if tid not in annoEnt.trInfo:
                            annoEnt.trInfo[tid]=dict()
                        if "trAlt" not in annoEnt.trInfo[tid]:
                            annoEnt.trInfo[tid]["trAlt"]=None
                        if re.search(r'^[ACGTN]+$', str(unify_a)) and str(rtidDetail['strd']) == '-':
                            annoEnt.trInfo[tid]['trAlt'] = BedAnno.BedAnno.rev_comp(unify_a)
                        else:
                            annoEnt.trInfo[tid]['trAlt'] = unify_a

                        tmp_trInfo = annoEnt.trInfo[tid]

                        tmp_trInfo['geneId'] = rtidDetail['gid']
                        tmp_trInfo['geneSym'] = rtidDetail['gsym']
                        tmp_trInfo['strd'] = rtidDetail['strd']
                        tmp_trInfo['primaryTag'] = rtidDetail['pr']

                        # assign left rna positions
                        # if no assignment then skip the right position assignment
                        if annoEnt.cal_hgvs_pos(tid=tid, tidDetail=rtidDetail, offset=total_left_ofst,LR=1):  # LR=1 left mode
                            if re.match(r'PROM', rtidDetail['blka']):
                                #      use P0 for sort
                                if 'trRefComp' not in tmp_trInfo:
                                    tmp_trInfo['trRefComp'] = dict()
                                if 'P0' not in tmp_trInfo['trRefComp']:
                                    tmp_trInfo['trRefComp']['P0']=list()
                                tmp_trInfo['trRefComp']['P0'].insert(0, 0)
                            elif re.match(r'IVS', rtidDetail['exin']):
                                if 'trRefComp' not in tmp_trInfo:
                                    tmp_trInfo['trRefComp']=dict()
                                if rtidDetail['exin'] not in tmp_trInfo['trRefComp']:
                                    tmp_trInfo['trRefComp'][rtidDetail['exin']]=list()
                                tmp_trInfo['trRefComp'][rtidDetail['exin']].insert(0, 0)
                            elif not strd:
                                m = re.match(r'\+(\d+)', str(tmp_trInfo['rnaEnd']))
                                if m:
                                    if 'trRefComp' not in tmp_trInfo:
                                        tmp_trInfo['trRefComp']=dict()
                                    if 'Z999' not in tmp_trInfo['trRefComp']:
                                        tmp_trInfo['trRefComp']['Z999']=list()
                                    tmp_trInfo['trRefComp']['Z999'].append(0)
                                    tmp_trInfo['trRefComp']['Z999'].append(m.group(1))

                            # let right rna positions to record Exon block

                            if total_left_ofst > 0:
                                rpreLeft[tid] = annoEnt.cal_hgvs_pos(tid=tid, tidDetail=rtidDetail, offset=total_left_ofst,LR=0, noassign=1)  # LR=0  preLeft mode
                            elif tid not in rpreLeft:

                                # no left block of annotation for this tid
                                # then 5'promoter left or 3'downstream
                                rpreLeft[tid] = 0

                            if not (tid in rpreLeft and type(rpreLeft[tid]) is dict):
                                # no left block of annotation for this tid
                                # then 5'promoter left or 3'downstream
                                if total_left_ofst == 0:
                                    preLeft_cDot = ""
                                    preLeft_nDot = ""
                                    preLeft_exin = ""
                                    preLeft_r = ""
                                    if strd:
                                        preLeft_exin = "."
                                        preLeft_r = "PROM"
                                        if tmp_trInfo['cdsBegin'] != "":
                                            m = re.match(r'(.*)-u(\d+)$', tmp_trInfo['cdsBegin'])
                                            if m:
                                                preLeft_cDot = "{}{}{}".format(m.group(1), "-u", (int(m.group(2)) + 1))
                                            else:  # for some sub region anno db entries
                                                preLeft_cDot = "{}{}".format("[uncertain preLeft]:",
                                                                             tmp_trInfo['cdsBegin'])
                                        m = re.match(r'\-?\d+$', tmp_trInfo['rnaBegin'])
                                        if tmp_trInfo['rnaBegin'] == 1:
                                            preLeft_nDot = -1
                                        elif m:
                                            preLeft_nDot = int(tmp_trInfo['rnaBegin']) - 1
                                        else:  # for some sub region anno db entries
                                            preLeft_nDot = "{}{}".format("[uncertain preLeft]:", tmp_trInfo['rnaBegin'])
                                    else:
                                        preLeft_exin = "."
                                        preLeft_r = "3D"
                                        preLeft_nDot = '+1'
                                        preLeft_cDot = '+d1'
                                    rpreLeft[tid] = {"nDot": preLeft_nDot, "cDot": preLeft_cDot, "exin": "preLeft_exin",
                                                     "r": preLeft_r}

                            if tid in rpreLeft and type(rpreLeft[tid]) is dict:
                                if strd:
                                    tmp_trInfo['preStart'] = rpreLeft[tid]
                                else:
                                    tmp_trInfo['postEnd'] = rpreLeft[tid]
                        else:
                            rpreLeft[tid] = annoEnt.cal_hgvs_pos(tid=tid, tidDetail=rtidDetail, offset=total_left_ofst,LR=0, noassign=1)  # LR=0 preLeft mode

                            continue

                    # $tid have been added into trInfo in the left end assignment
                    trinfoEnt = annoEnt.trInfo[tid]

                    if int(var.pos) == int(var.end):
                        # assume all insertion are on transcript
                        trinfoEnt['staInTr'] = 1
                        trinfoEnt['stoInTr'] = 1

                    # check if hit on transcript region
                    if int(rannodb[k]["sta"]) < int(rannodb[k]["sto"]) and int(var.pos) >= int(rannodb[k]["sta"]) and int(var.pos) < int(rannodb[k]["sto"]):
                        trinfoEnt['staInTr'] = 1
                    if int(rannodb[k]["sta"]) < int(rannodb[k]["sto"]) and int(var.end) > int(rannodb[k]["sta"]) and int(var.end) <= int(rannodb[k]["sto"]):
                        trinfoEnt['stoInTr'] = 1

                    # assign right rna positions
                    if (int(total_left_ofst) == 0) and (int(total_left_ofst) == int(total_right_ofst)) and (int(rtidDetail['wlen']) > 0) and (str(rtidDetail['mismatch']) == "") and (tid not in hit_badaln_ins):
                        # here is the part of insertion in edge,
                        # we don't need to recall cal_hgvs_pos,
                        # just use the already generated infomation
                        # About genepart judgement for this edge insertion case
                        # 1. promoter-any     => promoter
                        # 2. any-3'downstream => 3'downstream
                        # 3. utr-cds          => utr
                        # 4. exon-intron      => exon
                        if strd:
                            trinfoEnt['rnaEnd'] = rpreLeft[tid]["nDot"]
                            trinfoEnt['cdsEnd'] = rpreLeft[tid]["cDot"]
                            trinfoEnt['ei_End'] = rpreLeft[tid]["exin"]
                            trinfoEnt['r_End'] = rpreLeft[tid]["r"]
                        else:
                            trinfoEnt['rnaBegin'] = rpreLeft[tid]["nDot"]
                            trinfoEnt['cdsBegin'] = rpreLeft[tid]["cDot"]
                            trinfoEnt['ei_Begin'] = rpreLeft[tid]["exin"]
                            trinfoEnt['r_Begin'] = rpreLeft[tid]["r"]

                        # correct genepartSO for different type of edge insertion
                        if trinfoEnt['r_Begin'] == 'PROM' or trinfoEnt['r_End'] == 'PROM':
                            trinfoEnt['genepartSO'] = '167'
                            trinfoEnt['trRefComp']['P0'] = [0, 0]
                            trinfoEnt['exin'] = '.'
                            trinfoEnt['r'] = 'PROM'
                        elif trinfoEnt['r_Begin'] == '3D' or trinfoEnt['r_End'] == '3D':
                            # 3'downstream insertion will be ignored as intergenic variantion
                            del annoEnt.trInfo[tid]
                            continue
                        elif (re.search(r"EX", trinfoEnt['ei_End']) and re.search(r"IVS", trinfoEnt['ei_Begin'])) or (
                            re.search(r"^[53]U", trinfoEnt['r_End']) and re.search(r"^C", trinfoEnt['r_Begin'])):
                            trinfoEnt['genepartSO'] = BedAnnoAnno.getSOfromR(trinfoEnt['r_End'])
                            if "trRefComp" not in trinfoEnt or trinfoEnt['ei_End'] not in trinfoEnt['trRefComp']:
                                if "trRefComp" not in trinfoEnt:
                                    trinfoEnt['trRefComp']=dict()
                                if trinfoEnt['ei_End'] not in trinfoEnt['trRefComp']:
                                    trinfoEnt['trRefComp'][trinfoEnt['ei_End']]=None
                                trinfoEnt['trRefComp'][trinfoEnt['ei_End']] = 0
                            trinfoEnt['exin'] = trinfoEnt['ei_End']
                            trinfoEnt['r'] = trinfoEnt['r_End']
                        elif (re.search(r"EX", trinfoEnt['ei_Begin']) and re.search(r"IVS", trinfoEnt['ei_End'])) or (
                            re.search(r"^[53]U", trinfoEnt['r_Begin']) and re.search(r"^C", trinfoEnt['r_End'])):
                            trinfoEnt['genepartSO'] = BedAnnoAnno.getSOfromR(trinfoEnt['r_Begin'])
                            if "trRefComp" not in trinfoEnt or trinfoEnt['ei_End'] not in trinfoEnt['trRefComp']:
                                if "trRefComp" not in trinfoEnt:
                                    trinfoEnt['trRefComp'] = dict()
                                if trinfoEnt['ei_Begin'] not in trinfoEnt['trRefComp']:
                                    trinfoEnt['trRefComp'][trinfoEnt['ei_Begin']]=None
                                trinfoEnt['trRefComp'][trinfoEnt['ei_Begin']] = 0
                            trinfoEnt['exin'] = trinfoEnt['ei_Begin']
                            trinfoEnt['r'] = trinfoEnt['r_Begin']
                        else:  # no genepart change needed
                            if trinfoEnt['ei_Begin'] == trinfoEnt['ei_End']:
                                if re.match(r'EX', trinfoEnt['ei_Begin']):
                                    if "trRefComp" not in trinfoEnt or trinfoEnt['ei_Begin'] not in trinfoEnt[
                                        'trRefComp']:
                                        trinfoEnt['trRefComp'][trinfoEnt['ei_Begin']] = 0
                                elif trinfoEnt['ei_Begin'] == '?':
                                    trinfoEnt['trRefComp']['Q1'] = 0
                                else:
                                    trinfoEnt['trRefComp'][trinfoEnt['ei_Begin']] = [0, 0]
                            trinfoEnt['exin'] = trinfoEnt['ei_Begin']
                            trinfoEnt['r'] = trinfoEnt['r_Begin']
                    elif annoEnt.cal_hgvs_pos(tid=tid, tidDetail=rtidDetail, offset=total_right_ofst,LR=0):  # LR=0 right mode
                        # assign the trRefComp
                        tmp_exin = rtidDetail['exin']
                        if re.match(r'EX', tmp_exin):
                            if (tid not in tidExblk) or (tmp_exin not in tidExblk[tid]) or ("{}{}{}".format(rtidDetail['nsta'], ',', rtidDetail['nsto']) not in tidExblk[tid][tmp_exin]):
                                # non Insertion on refseq, and first occured exon
                                if tid not in tidExblk:
                                    tidExblk[tid]=dict()
                                if tmp_exin not in tidExblk[tid]:
                                    tidExblk[tid][tmp_exin]=dict()
                                tidExblk[tid][tmp_exin][
                                    "{}{}{}".format(rtidDetail['nsta'], ',', rtidDetail['nsto'])] = 1
                                cmpStart, cmpEnd = None, None
                                if strd:
                                    # select the larger one
                                    if 0 < BedAnno.BedAnno.cmpPos(trinfoEnt['rnaBegin'], rtidDetail['nsta']):
                                        cmpStart = rtidDetail['nsta']
                                    else:
                                        cmpStart = trinfoEnt['rnaBegin']
                                    # select the smaller one
                                    if 0 < BedAnno.BedAnno.cmpPos(trinfoEnt['rnaEnd'], rtidDetail['nsto']):
                                        cmpEnd = trinfoEnt['rnaEnd']
                                    else:
                                        cmpEnd = rtidDetail['nsto']
                                else:
                                    # select the larger one
                                    if 0 < BedAnno.BedAnno.cmpPos(trinfoEnt['rnaBegin'], rtidDetail['nsto']):
                                        cmpStart = rtidDetail['nsto']
                                    else:
                                        cmpStart = trinfoEnt['rnaBegin']
                                    # select the smaller one
                                    if 0 < BedAnno.BedAnno.cmpPos(trinfoEnt['rnaEnd'], rtidDetail['nsta']):
                                        cmpEnd = trinfoEnt['rnaEnd']
                                    else:
                                        cmpEnd = rtidDetail['nsta']
                                if (not re.match(r'\d+$', str(cmpStart))) or (not re.match(r'\d+$', str(cmpEnd))):
                                    raise Error.Error("{}{}{}{}{}".format("Error: unknown bug cmpStart [",cmpStart , "], cmpEnd [", cmpEnd, "]"))
                                if strd and "trRefComp" in trinfoEnt and "Z999" in trinfoEnt['trRefComp']:
                                    del trinfoEnt['trRefComp']['Z999']
                                if 'trRefComp' not in trinfoEnt:
                                    trinfoEnt['trRefComp']=dict()
                                if tmp_exin not in trinfoEnt['trRefComp']:
                                    trinfoEnt['trRefComp'][tmp_exin]=0
                                trinfoEnt['trRefComp'][tmp_exin] = int(trinfoEnt['trRefComp'][tmp_exin]) + int(cmpEnd) - int(cmpStart) + 1
                                m = re.match(r'\+(\d+)', str(trinfoEnt['rnaEnd']))
                                if m and strd and re.search(r'E$', str(tmp_exin)):
                                    if 'Z999' not in trinfoEnt['trRefComp']:
                                        trinfoEnt['trRefComp']['Z999']=list()
                                    trinfoEnt['trRefComp']['Z999'].insert(0 , int(unify_rl) - int(m.group(1)))
                                    trinfoEnt['trRefComp']['Z999'].insert(1 ,int(unify_rl))
                        elif rtidDetail['blka'] == '?':
                            # problem for trRef assign for new version
                            af_sta, af_sto = None, None
                            if int(rannodb[k]['sta']) < int(unify_p):
                                af_sta = unify_p
                            else:
                                af_sta = rannodb[k]['sta']
                            if int(rannodb[k]['sto']) < int(unify_p) + int(unify_rl):
                                af_sto = rannodb[k]['sto']
                            else:
                                af_sto = int(unify_p) + int(unify_rl)
                            trinfoEnt['trRefComp']['Q1'] = int(trinfoEnt['trRefComp']['Q1']) + (int(af_sto) - int(af_sta))
                        else:  # for intron region and PROM
                            if re.match(r'PROM', rtidDetail['blka']):
                                tmp_exin = 'P0'  # for sort
                            if strd and "trRefComp" in trinfoEnt and "Z999" in trinfoEnt['trRefComp']:
                                del trinfoEnt['trRefComp']["Z999"]
                            if "trRefComp" not in trinfoEnt or tmp_exin not in trinfoEnt['trRefComp']:
                                if "trRefComp" not in trinfoEnt:
                                    trinfoEnt['trRefComp']=dict()
                                if tmp_exin not in trinfoEnt['trRefComp']:
                                    trinfoEnt['trRefComp'][tmp_exin]=list()
                                trinfoEnt['trRefComp'][tmp_exin].insert(0 ,int(rannodb[k]['sta']) - int(unify_p))
                            less = None
                            if int(rannodb[k]['sto']) > (int(unify_p) + int(unify_rl)):
                                less = int(unify_p) + int(unify_rl)
                            else:
                                less = rannodb[k]['sto']
                            trinfoEnt['trRefComp'][tmp_exin].insert(1,int(less) - int(unify_p))
                            if str(tmp_exin) == 'P0' and not strd:
                                trinfoEnt['trRefComp'][tmp_exin][1] = unify_rl

                    rpostRight = annoEnt.cal_hgvs_pos(tid=tid, tidDetail=rtidDetail, offset=total_right_ofst, LR=1,noassign=1)  # LR=1 post right mode
                    if type(rpostRight) is dict:
                        if strd:
                            trinfoEnt['postEnd'] = rpostRight
                        else:
                            trinfoEnt['preStart'] = rpostRight
                            # debug
                            #                print STDERR "AnnoEnt Check: "
                            #                  . Dumper( $annoEnt->{trInfo}->{$tid} )
                            #                  if (  exists $annoEnt->{trInfo}
                            #                    and exists $annoEnt->{trInfo}->{$tid} );

        # debug
        #    print STDERR "Annotations before check and uniform:", Dumper($annoEnt);

        # final check and uniform trinfo
        if hasattr(annoEnt, "trInfo"):
            for t in annoEnt.trInfo.keys():
                tinfo = annoEnt.trInfo[t]
                if "staInTr" in tinfo and "stoInTr" in tinfo:
                    del tinfo['staInTr']
                    del tinfo['stoInTr']
                elif "staInTr" not in tinfo and "stoInTr" not in tinfo and "preStart" in tinfo and tinfo['preStart'][
                    'r'] != 'PROM':
                    # hit 3'downstream only
                    del annoEnt.trInfo[t]
                    continue
                else:
                    if "staInTr" in tinfo:
                        del tinfo['staInTr']
                    if "stoInTr" in tinfo:
                        del tinfo['stoInTr']
                if (re.search(r'\+', tinfo['strd']) and ("rnaEnd" not in tinfo)) or (re.search(r'\-', tinfo['strd']) and "rnaBegin" not in tinfo):
                    del annoEnt.trInfo[t]
                elif "cdsEnd" in tinfo and tinfo['cdsEnd'] == '':
                    # Downstream annotation fail
                    tinfo['cdsBegin'] = ''
            if 0 == len(annoEnt.trInfo.keys()):
                del annoEnt.trInfo

        # debug
        #    print STDERR "DEBUG: Done getTrPosition\n";

        if not new_aeIndex:
            new_aeIndex = aeIndex

        return new_aeIndex