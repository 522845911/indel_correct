#!/usr/bin/env python
# coding utf-8

import re, copy, sys, json,BedAnno,Error
sys.dont_write_bytecode = True


class BedAnnoVar(object):
    MAX_COMPLEX_PARSING = 200


    # =head1 METHOD
    #
    # =head2 new
    #
    # 	About   : Create a new object class entry, BedAnno::Var,
    # 			  parse the variation directly by the ref and alt string.
    # 	Usage   : my $var = BedAnnoVar( chr, start, end, ref, alt );
    # 		   or my $var = BedAnnoVar( chr, pos, ref, alt );
    # 		   or my $var = BedAnnoVar( varInput_dict );
    # 	Args    : Input can be variable format
    # 			  1. 5 parameters format: CG shell list format: chr,start,end,ref,alt
    # 			  2. 4 parameters format: VCF format: chr,pos,ref,alt
    # 			  3. Crawler input object: A dict with nessesary keys:
    # 				 chr,begin,referenceSequence,variantSequence,
    # 				 optional key is "end", if end specified,
    # 				 coordinates are treat as 0-based, otherwise, use 1-based (VCF)
    # 	Returns : a new BedAnnoVar instance :
    # 			{
    # 				chr    => chr,
    # 				pos    => start,          # 0-based start
    # 				end    => end,
    # 				ref    => ref,
    # 				alt    => alt,
    # 				reflen => ref_len,
    # 				altlen => alt_len,        # not exists if no-call
    # 				guess  => varType,        # the output varType
    # 				imp    => imp_varType,    # the implicit varType
    # 				sm     => sm,             # single/multiple base indicator
    # 										   # equal/non-equal length indicator
    #
    # 				# for hgvs naming convinient, reparse delins(guess),
    # 				# in forward and reverse strand separately,
    # 				# If the result are the same, then only
    # 				# give the following optional
    # 				# rescaled strand-same description group
    # 				bp  => bc_pos,       # backward compatible pos, 0-based
    # 				br  => bc_ref,       # backward compatible ref string
    # 				ba  => bc_alt,       # backward compatible alt string
    # 				brl => bc_reflen,    # backward compatible ref length
    # 				bal => bc_altlen,    # backward compatible alt length
    #
    # 				# otherwise, the following '+', '-',
    # 				# structure will be generated to reflect
    # 				# the difference. they are all optional
    #
    # 				'+' => {
    #
    # 				  # This group simplely trim off the leading same chars
    # 				  # on forward strand, and then trim the same tail
    # 				  bp  => backward_fpos,
    # 				  br  => backward_fref,
    # 				  ba  => backward_falt,
    # 				  brl => backward_freflen,
    # 				  bal => backward_faltlen,
    #
    # 				},
    #
    # 				'-' => {
    # 				  # similar to '+', but for reverse strand
    # 				},
    #
    # 				# this group gives ref/alt string based on the rule
    # 				# with 'rep' annotation available
    # 				p      => rep_left_pos,         # repeat pos, 0-based
    # 				r      => rep_ref,              # repeat ref string
    # 				a      => rep_alt,              # repeat alt string
    # 				rl     => rep_reflen,           # repeat ref length
    # 				al     => rep_altlen,           # repeat alt length
    # 				rep    => repeat_element,
    # 				replen => repeat_element_length,
    # 				ref_cn => copy_number_in_ref,
    # 				alt_cn => copy_number_in_alt,
    #
    # 				# for equal length long substitution
    # 				# record the separated snvs positions
    # 				# all positions are 1 based.
    # 				sep_snvs => [ snv_pos1, snv_pos2, ... ],
    # 			}
    #
    # =cut
    def __init__(self, *args):  # chr, start, end, ref, alt
        chr, start, end, ref, alt = None, None, None, None, None
        var = dict()

        if (len(args) == 1 and isinstance(args[0], (dict))):
            var_dict = copy.deepcopy(args[0])
            if (("chr" not in var_dict) or ("begin" not in var_dict) or (
                    "end" not in var_dict and "referenceSequence" not in var_dict)):
                raise Error.Error("Error: unavailable object. need keys: " + "chr, start, alt, ref specified.")

            var = copy.deepcopy(var_dict)

            chr = var_dict["chr"]
            start = var_dict["begin"]
            if "end" in var_dict:
                end = var_dict["end"]
            if ("variantSequence" in var_dict and var_dict["variantSequence"] is not None):
                alt = var_dict["variantSequence"]
            else:
                alt = ""
            if (re.search("^null$", alt, re.IGNORECASE)):
                alt = ""

            if ("referenceSequence" in var_dict and var_dict["referenceSequence"] is None):
                var_dict["referenceSequence"] = ""

            if ("end" in var_dict and start == end):
                ref = ""
            else:
                if ("referenceSequence" in var_dict):
                    ref = var_dict["referenceSequence"]
                else:
                    ref = "="
            if (re.search("^null$", ref, re.IGNORECASE)):
                ref = ""

            # clean var dict
            del var["chr"]
            del var["begin"]

            if ("end" in var):
                del var["end"]

            if ("referenceSequence" in var):
                del var["referenceSequence"]

            if ("variantSequence" in var):
                del var["variantSequence"]
        else:
            if (len(args) < 4):
                raise Error.Error("Error: not enough args, need at least 4 args.")
            if (len(args) == 4):
                chr, start, end, ref = args
            elif (len(args) == 5):
                chr, start, end, ref, alt = args

        if end is None or (not re.match("^\d+$", str(end))):  # from VCF v4.1 1-based start
            if (end is not None):
                alt = ref
                ref = end
            try:
                ref = BedAnnoVar.normalise_seq(ref)
                alt = BedAnnoVar.normalise_seq(alt)
            except Exception,e:
                raise e
            #if ((ref != alt or len(ref) > 1) and (ref[0] == alt[0])):
            if ref != alt and len(ref) >=1 and len(alt)>=1 and ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
            else:
                start =int(start)- 1  # change to 0-based start
            rl = len(ref)
            end = int(start) + rl

        len_ref = int(end) - int(start)  # chance to annotate long range
        varType, implicit_varType, sm = BedAnnoVar.guess_type(len_ref, ref, alt)

        chr = BedAnnoVar.normalise_chr(chr)
        update_dict = {"chr": chr, "pos": start, "ref": ref.upper(), "end": end, "alt": alt.upper(), "reflen": len_ref,
                       "guess": varType, "imp": implicit_varType, "sm": sm}
        var.update(update_dict)  # remain some extra keys in the var hash, e.g. var_id etc.

        if alt != '?':
            var["altlen"] = len(alt)

        for i in var:
            setattr(self, i, var[i])

        if not ((varType == 'no-call' or implicit_varType != 'delins') or (
                len_ref > BedAnnoVar.MAX_COMPLEX_PARSING or len(alt) > BedAnnoVar.MAX_COMPLEX_PARSING)):
            self.parse_complex()

    @staticmethod
    def normalise_seq(seq):
        seq, _ = re.subn(r"\s+", "", seq)
        if seq == '.':
            seq = ""
        seq = seq.upper()
        m = re.search(r"[^ACGTN]", seq)
        if m:
            raise Error.Error("Error: unrecognized pattern exists,no multiple alts please. [" + seq + "]")
        return seq

    @staticmethod
    def TO_JSON(data):
        return json.dumps(data)

    @staticmethod
    def normalise_chr(chr):
        chr = str(chr)
        chr, _ = re.subn("^chr", "", chr, re.IGNORECASE)
        if (chr == '23'):
            chr = 'X'
        elif (chr == '24'):
            chr = 'Y'
        elif (chr == '25'):
            chr = 'MT'
        elif (re.search(r"^[xym]$", chr, re.IGNORECASE)):
            chr = chr.upper()
            if chr == 'M':
                chr += 'T'
        elif (re.search(r"^M_", chr, re.IGNORECASE)):
            chr = 'MT'
        return chr

    # guess only by length for get_internal result
    # the only delins without no-call left
    # get_internal will recognize the real mutant
    # region, discard the influnce by other
    # sample's result.
    @staticmethod
    def guess_type_by_length(rl, al):
        imp_guess, sm = "", ""
        if (rl == 0):
            imp_guess = 'ins'
            sm = 0
        elif (rl == 1):
            if (al == 1):
                imp_guess = 'snv'
            elif (al == 0):
                imp_guess = 'del'
            else:
                imp_guess = 'delins'
            sm = 1
        else:
            imp_guess = 'del' if al == 0 else "delins"
            sm = 3 if rl == al else 2

        return imp_guess, sm

    # imp_varType: implicit variant type is for HGVS naming
    # varType    : to be output as varType name. can be as the key
    #              to get the SO term by Name2SO hash.
    # sm         : single or multiple bases tag
    #              0 - for insertion, 0 base
    #              1 - for single base variants
    #              2 - for non-equal-length multiple bases variants
    #              3 - for equal-length multiple bases delins
    @staticmethod
    def guess_type(len_ref, ref, alt):
        imp_varType, varType, sm = "", "", ""
        if (len_ref == 0):
            imp_varType = 'ref' if ref == alt else 'ins'
            sm = 0
        elif (len_ref == 1):
            if (ref == alt):
                imp_varType = 'ref'
            elif (alt == ''):
                imp_varType = 'del'
            elif (1 == len(alt)):
                imp_varType = 'snv'
            else:
                imp_varType = 'delins'
            sm = 1
        elif (len_ref > 1):
            if (ref == alt):
                imp_varType = 'ref'
            elif (alt == ''):
                imp_varType = 'del'
            else:
                imp_varType = 'delins'
            sm = 2 if (len(ref) != len(alt)) else 3  # sm = 2 non-equal-length delins
            # sm = 3     # equal-length subs
        varType = 'no-call' if (alt == '?' or alt == 'N') else imp_varType
        return varType, imp_varType, sm

    # check whether the smaller with inserted $cn copies repeat elements
    # is the same with larger one
    @staticmethod
    def check_insrep(larger, smaller, repeat, cn):
        ind = smaller.find(repeat)
        if (ind > -1):
            temp = smaller
            temp = temp[:ind] + repeat * cn + temp[ind:]
            if (larger == temp):
                return 1
        return 0

    # =head2 get_internal
    #
    # 	About   : recalculate ref alt for delins and mutiple sample caused ref-alt pair.
    # 			  depend on different strand of the to-annotated gene or transcript,
    # 			  the offset may be different for the same ref and alt,
    # 			  because of the 3'end nearest annotation rules.
    # 	Usage   : my $rephase = get_internal( $ref, $reflen, $alt, $altlen );
    # 	Returns : a hash ref of :
    # 				{
    # 					'+' => $f_lofs,
    # 					'-' => $r_lofs,
    # 					'r' => $new_ref_len,
    # 					'a' => $new_alt_len
    # 				}
    #
    # =cut
    @staticmethod
    def get_internal(ref, reflen, alt, altlen):
        shorter = reflen if reflen < altlen else altlen
        lgo, loff, rgo, roff = 1, 0, 1, 0
        for i in range(0, shorter):
            if (lgo and ref[i] == alt[i]):
                loff += 1
            else:
                lgo = 0
            if (rgo and ref[-(i + 1)] == alt[-(i + 1)]):
                roff += 1
            else:
                rgo = 0
            if lgo == 0 and rgo == 0:
                break
        new_ref_len, new_alt_len = None, None
        if (shorter >= loff + roff):
            new_ref_len = reflen - loff - roff
            new_alt_len = altlen - loff - roff
            return {'+': loff, '-': loff, 'r': new_ref_len, 'a': new_alt_len}
        else:
            new_ref_len = reflen - shorter
            new_alt_len = altlen - shorter
            return {'+': loff, '-': (shorter - roff), 'r': new_ref_len, 'a': new_alt_len}

    # check sign for diff-array
    @staticmethod
    def check_sign(rc):
        if ((rc[0] > 0 and len([i for i in rc if i < 0]) > 0) or (rc[0] < 0 and len([i for i in rc if i > 0]) > 0)):
            return 0
        return 1

    @staticmethod
    def count_content(s):
        s = s.upper()
        l = len(s)
        count = [l] + [0] * BedAnno.BedAnno.AAcount
        for char in s:
            if (char not in BedAnno.BedAnno.AAnumber):
                raise Error.Error("unknown base [" + char + "]")
            count[BedAnno.BedAnno.AAnumber[char]] += 1
        return count

    # check length and content consistent-divisability for string and diff-array
    @staticmethod
    def check_div(s, rc):
        rcs = BedAnnoVar.count_content(s)
        if rc[0] % rcs[0] != 0:
            return 0
        div = rc[0] / rcs[0]
        for i in range(1, BedAnno.BedAnno.AAcount + 1):
            if (rc[i] != rcs[i] * div):
                return 0
        return div

    # =head2 getUnifiedVar
    #
    # 	About   : uniform the pos and ref/alt pair selection,
    # 			  after new(), give info for HGVS naming.
    # 	Usage   : my @unified_desc = $var->getUnifiedVar($strd);
    # 	Args    : BedAnno::Var entry and current strand for annotation.
    # 	Returns : an array of (
    # 				$pos,  # 0-based start pos
    # 				$ref,  # reference bases
    # 				$alt,  # called bases
    # 				$reflen, # reference len
    # 				$altlen )# called len, undef if no-call
    #
    # =cut
    def getUnifiedVar(self, strd):
        consPos = self.pos
        consRef = self.ref
        consAlt = self.alt
        consRL = self.reflen

        if ("altlen" not in self.__dict__):
            return (consPos, consRef, consAlt, consRL,None)

        consAL = self.altlen

        if ("p" in self.__dict__):  # rep
            consPos = self.p
            consRef = self.r
            consAlt = self.a
            consRL = self.rl
            consAL = self.al
        elif ("bp" in self.__dict__):  # complex bc strand same
            consPos = self.bp
            consRef = self.br
            consAlt = self.ba
            consRL = self.brl
            consAL = self.bal
        elif ("strd" in self.__dict__):  # complex bc strand diff
            consPos = self.strd["bp"]
            consRef = self.strd["br"]
            consAlt = self.strd["ba"]
            consRL = self.strd["brl"]
            consAL = self.strd["bal"]

        return (consPos, consRef, consAlt, consRL, consAL)

    def get_gHGVS(self):
        gHGVS = 'g.'
        if (re.match("M", self.chr)):  # hit mito chromosome
            gHGVS = 'm.'

        imp = self.imp
        sm = self.sm
        (pos, ref, alt, reflen, altlen) = self.getUnifiedVar('+')
        pos, ref, alt, reflen = str(pos), str(ref), str(alt), str(reflen)
        if (imp == 'snv'):
            gHGVS = "{}{}{}{}{}".format(gHGVS, pos, ref, '>', alt)
        elif (imp == 'ins'):
            gHGVS = "{}{}{}{}{}{}".format(gHGVS, pos, '_', (int(pos) + 1), 'ins', alt)
        elif (str(imp) == 'del' or str(imp) == 'delins'):
            gHGVS = "{}{}".format(gHGVS, (int(pos) + 1))
            if (int(sm) > 1):
                gHGVS = "{}{}{}".format(gHGVS, '_', (int(pos) + int(reflen)))
            if not re.search("^[ACGTN]+$", ref):
                ref = ""
            gHGVS = "{}{}{}".format(gHGVS, 'del', ref)
            if ("delins" == imp):
                gHGVS += ('ins' + alt)
        elif (str(imp) == 'rep'):
            # 1bp del/delins
            gHGVS = "{}{}".format(gHGVS, (int(pos) + 1))
            if (int(self.ref_cn) == 1 and int(self.alt_cn) == 2):  # dup
                if (int(sm) > 1):
                    gHGVS = "{}{}{}".format(gHGVS, '_', int(pos) + int(self.replen))
                gHGVS = "{}{}{}".format(gHGVS, 'dup', self.rep)
            else:
                gHGVS = "{}{}{}{}{}{}{}".format(gHGVS, self.rep, '[', self.ref_cn, '>', self.alt_cn, ']')
        elif (str(imp) == 'ref'):
            if (int(reflen) == 1):
                gHGVS = "{}{}{}{}".format(gHGVS, pos, ref, '=')
            else:
                gHGVS = "{}{}{}{}{}".format(gHGVS, pos, '_', int(pos) + int(reflen) - 1, '=')
        else:
            raise Error.Error("Can not recognize type " + str(imp) + ".")
        self.gHGVS = gHGVS

    def parse_complex(self):
        ref, alt, len_ref, len_alt = self.ref, self.alt, self.reflen, self.altlen
        get_rst = self.get_internal(ref, len_ref, alt, len_alt)
        if (get_rst["r"] != len_ref):
            imp_guess, sm = BedAnnoVar.guess_type_by_length(get_rst["r"], get_rst["a"])
            self.imp = imp_guess
            self.sm = sm
            if (get_rst["+"] != get_rst["-"]):
                for strd in ("+", "-"):
                    setattr(self, strd, dict())
                    getattr(self, strd)["bp"] = self.pos + get_rst[strd]
                    getattr(self, strd)["brl"] = get_rst["r"]
                    getattr(self, strd)["bal"] = get_rst["a"]
                    getattr(self, strd)["br"] = self.ref[get_rst[strd]:get_rst[strd] + get_rst["r"]]
                    getattr(self, strd)["ba"] = self.alt[get_rst[strd]:get_rst[strd] + get_rst["a"]]
            else:
                self.bp = self.pos + get_rst["+"]
                self.brl = get_rst["r"]
                self.bal = get_rst["a"]
                self.br = self.ref[get_rst["+"]:get_rst["+"] + get_rst["r"]]
                self.ba = self.alt[get_rst["+"]:get_rst["+"] + get_rst["a"]]

        if (self.sm == 3):  # equal length long subs
            separate_snvs = list()
            for i in range(0, self.reflen):
                if (self.ref[i] == self.alt[i]):
                    continue
                separate_snvs.append(self.pos + i + 1)  # 1-based
            self.sep_snvs = separate_snvs

        rc_ref = BedAnnoVar.count_content(ref)
        rc_alt = BedAnnoVar.count_content(alt)
        diff = list()
        for i in range(0, BedAnno.BedAnno.AAcount + 1):
            diff.append(rc_ref[i] - rc_alt[i])

        # check if the sign of all base diff are consistent.
        if (BedAnnoVar.check_sign(diff)):  # possible short tandom repeat variation
            absdiff = map(abs, diff)
            larger, smaller, llen, slen = None, None, None, None
            if (len_ref > len_alt):
                larger = ref
                llen = len_ref
                smaller = alt
                slen = len_alt
            else:
                larger = alt
                llen = len_alt
                smaller = ref
                slen = len_ref

            has = dict()
            for rep in range(llen, 0, -1):
                for m in re.finditer(r"([A-Z]+)(?:\1){" + str(rep) + "}", larger):
                    if m.group(1) in has:
                        continue
                    rep_el = m.group(1)
                    lofs = m.start()  # the prematched string length

                    cn = BedAnnoVar.check_div(rep_el, absdiff)
                    if (cn and BedAnnoVar.check_insrep(larger, smaller, rep_el, cn)):
                        lenrep = len(rep_el)

                        self.p, self.rep, self.replen = self.pos + lofs, rep_el, lenrep
                        rep += 1  # add the first copy of element

                        l_cn = rep
                        s_cn = rep - cn
                        l_wholerep = str(rep_el) * l_cn
                        s_wholerep = str(rep_el) * s_cn
                        l_replen = int(lenrep) * l_cn
                        s_replen = int(lenrep) * s_cn

                        if (llen == len_ref):  # ref is longer
                            self.ref_cn, self.alt_cn, self.r, self.a, self.rl, self.al = l_cn, s_cn, l_wholerep, s_wholerep, l_replen, s_replen
                        else:  # alt is longer
                            self.alt_cn, self.ref_cn, self.a, self.r, self.al, self.rl = l_cn, s_cn, l_wholerep, s_wholerep, l_replen, s_replen

                        self.imp = 'rep'
                        self.sm = 1 if self.rl == 1 else 2
                        return
                    has[rep_el] = 1
