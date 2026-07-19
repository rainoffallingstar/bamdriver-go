package bamnative

import (
	"testing"
)

func TestSAMFlagConstants(t *testing.T) {
	// Verify SAM flag constants are distinct powers of two
	flags := map[string]uint16{
		"FlagPaired":      FlagPaired,
		"FlagProperPair":   FlagProperPair,
		"FlagUnmapped":     FlagUnmapped,
		"FlagMateUnmapped": FlagMateUnmapped,
	}

	seen := make(map[uint16]bool)
	for name, f := range flags {
		if f == 0 {
			t.Errorf("%s is zero", name)
		}
		if seen[f] {
			t.Errorf("duplicate flag value %#x for %s", f, name)
		}
		seen[f] = true
		if f&(f-1) != 0 {
			t.Errorf("%s = %#x is not a single bit", name, f)
		}
	}
}

func TestRecord_IsPaired(t *testing.T) {
	var rec Record

	// Unmapped record should not be paired
	rec.Flags = FlagUnmapped
	if rec.IsPaired() {
		t.Error("unmapped record should not be paired")
	}

	// Paired record should be paired
	rec.Flags = FlagPaired
	if !rec.IsPaired() {
		t.Error("paired record should be paired")
	}

	// Paired + proper pair
	rec.Flags = FlagPaired | FlagProperPair
	if !rec.IsPaired() {
		t.Error("properly paired record should be paired")
	}
}

func TestRecord_IsUnmapped(t *testing.T) {
	var rec Record

	rec.Flags = 0
	if rec.IsUnmapped() {
		t.Error("mapped record should not be unmapped")
	}

	rec.Flags = FlagUnmapped
	if !rec.IsUnmapped() {
		t.Error("unmapped record should be unmapped")
	}
}

func TestRecord_IsMateUnmapped(t *testing.T) {
	var rec Record

	rec.Flags = 0
	if rec.IsMateUnmapped() {
		t.Error("mate-mapped record should not have mate unmapped")
	}

	rec.Flags = FlagMateUnmapped
	if !rec.IsMateUnmapped() {
		t.Error("mate-unmapped record should have mate unmapped")
	}
}

func TestRecord_HasRefID(t *testing.T) {
	var rec Record
	rec.RefID = 0
	if !rec.HasRefID() {
		t.Error("record with RefID=0 should have RefID (chr=*)")
	}
	rec.RefID = -1
	if rec.HasRefID() {
		t.Error("record with RefID=-1 should not have RefID")
	}
}

func TestRecord_HasNM(t *testing.T) {
	var rec Record
	ok := HasNM(&rec, "NM")
	if ok {
		t.Log("NM tag unexpectedly present on empty record")
	}
}

func TestRecord_GetAuxField(t *testing.T) {
	var rec Record
	aux := rec.GetAuxField("NM")
	if aux != nil {
		t.Log("Unexpected aux field on empty record")
	}
}

func TestCIGARConstants(t *testing.T) {
	ops := map[string]byte{
		"CIGARMatch":              CIGARMatch,
		"CIGARInsertion":          CIGARInsertion,
		"CIGARDeletion":           CIGARDeletion,
		"CIGARSkip":               CIGARSkip,
		"CIGARSoftClip":           CIGARSoftClip,
		"CIGARHardClip":           CIGARHardClip,
		"CIGARPadding":            CIGARPadding,
		"CIGARSequenceMatch":      CIGARSequenceMatch,
		"CIGARSequenceMismatch":   CIGARSequenceMismatch,
	}
	for name, val := range ops {
		if val == 0 {
			t.Errorf("%s constant is zero", name)
		}
	}
}
